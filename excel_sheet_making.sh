## Command for installing dependencies 
conda create -n ncbi-tools -c conda-forge -c bioconda entrez-direct
#
conda activate ncbi-tools
#
conda install -n ncbi-tools -c conda-forge perl
#
conda install -c conda-forge pandas openpyxl
#
conda activate ncbi-tools
#
which perl
perl -MTime::HiRes -e 'print "Perl is ready\n"'

## This is the script for the Excel sheet in just one go, and the main filters are Illumina, BGISEQ, and the WGS platform 

python3 - <<'PY'
import subprocess
import time
import xml.etree.ElementTree as ET
from shlex import quote
from urllib.parse import quote as urlquote

import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font
from openpyxl.utils import get_column_letter

ORGANISM = "Felis silvestris"

QUERY = (
    f'"{ORGANISM}"[Organism] AND '
    f'"strategy wgs"[Properties] AND '
    f'("platform illumina"[Properties] OR "platform bgiseq"[Properties])'
)

PREFIX = "felis_silvestris_illumina_bgiseq_wgs"
RUN_CSV = f"{PREFIX}_runs.csv"
BIOSAMPLE_XML = f"{PREFIX}_biosamples.xml"
OUTPUT_XLSX = f"{PREFIX}.xlsx"

BATCH_SIZE = 10
MAX_RETRIES = 4


def run(command):
    print(command)
    subprocess.run(
        ["bash", "-o", "pipefail", "-c", command],
        check=True,
    )


def capture(command):
    print(command)
    completed = subprocess.run(
        ["bash", "-o", "pipefail", "-c", command],
        check=True,
        text=True,
        stdout=subprocess.PIPE,
    )
    return completed.stdout


def normalise(name):
    return (name or "").strip().lower().replace(" ", "_").replace("-", "_")


def first_attribute(attributes, *names):
    for name in names:
        value = attributes.get(normalise(name), "")
        if value:
            return value
    return ""


def human_readable_count(value):
    if pd.isna(value):
        return ""

    try:
        value = int(value)
    except (TypeError, ValueError):
        return ""

    for divisor, unit in [
        (1_000_000_000, "G"),
        (1_000_000, "M"),
        (1_000, "K"),
    ]:
        if value >= divisor:
            return f"{value / divisor:.2f} {unit}"

    return f"{value:,}"


def fetch_biosample_batch(accessions):
    command = (
        f"efetch -db biosample "
        f"-id {quote(','.join(accessions))} "
        f"-format xml"
    )

    last_error = None

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            xml_text = capture(command)
            root = ET.fromstring(xml_text)

            if not root.findall(".//BioSample"):
                raise RuntimeError("No BioSample records returned in this batch.")

            return root

        except (
            subprocess.CalledProcessError,
            ET.ParseError,
            RuntimeError,
        ) as error:
            last_error = error

            if attempt == MAX_RETRIES:
                break

            wait_seconds = attempt * 3
            print(
                f"BioSample batch failed (attempt {attempt}/{MAX_RETRIES}). "
                f"Retrying in {wait_seconds} seconds..."
            )
            time.sleep(wait_seconds)

    raise RuntimeError(
        "NCBI could not return a valid BioSample batch after multiple attempts."
    ) from last_error


run(
    f"esearch -db sra -query {quote(QUERY)} "
    f"| efetch -format runinfo > {quote(RUN_CSV)}"
)

try:
    runs = pd.read_csv(RUN_CSV, dtype=str).fillna("")
except pd.errors.EmptyDataError as exc:
    raise RuntimeError(
        "No run metadata was returned. Check the organism name, "
        "Entrez Direct, and NCBI connectivity."
    ) from exc

if runs.empty:
    raise RuntimeError(
        "No matching WGS Illumina/BGISEQ runs were found for "
        f"{ORGANISM}."
    )

print(f"Run metadata rows retrieved: {len(runs)}")

biosample_ids = sorted({
    value.strip()
    for value in runs["BioSample"]
    if value and value.strip()
})

print(f"Unique BioSample accessions to retrieve: {len(biosample_ids)}")

combined_root = ET.Element("BioSampleSet")

for start in range(0, len(biosample_ids), BATCH_SIZE):
    batch = biosample_ids[start:start + BATCH_SIZE]

    print(
        f"Fetching BioSample batch "
        f"{start + 1}-{start + len(batch)} of {len(biosample_ids)}"
    )

    batch_root = fetch_biosample_batch(batch)

    for sample in batch_root.findall(".//BioSample"):
        combined_root.append(sample)

    time.sleep(0.4)

ET.ElementTree(combined_root).write(
    BIOSAMPLE_XML,
    encoding="utf-8",
    xml_declaration=True,
)

biosamples = {}

for sample in combined_root.findall(".//BioSample"):
    accession = sample.get("accession", "")
    title = sample.findtext("./Description/Title", default="").strip()
    attributes = {}

    for attribute in sample.findall("./Attributes/Attribute"):
        key = normalise(
            attribute.get("harmonized_name")
            or attribute.get("attribute_name")
        )
        value = (attribute.text or "").strip()

        if key and value:
            attributes.setdefault(key, value)

    biosamples[accession] = {
        "Description": title,
        "Breed": attributes.get("breed", ""),
        "Geographical location": first_attribute(
            attributes,
            "geo_loc_name",
            "geographic_location",
            "country",
        ),
        "Tissue": first_attribute(
            attributes,
            "tissue",
            "body_site",
            "isolation_source",
        ),
        "Sex": attributes.get("sex", ""),
        "Birth date": first_attribute(
            attributes,
            "birth_date",
            "birthdate",
            "date_of_birth",
            "date of birth",
            "dob",
        ),
        "Death date": first_attribute(
            attributes,
            "death_date",
            "deathdate",
            "date_of_death",
            "date of death",
            "dod",
        ),
        "Disease": first_attribute(
            attributes,
            "disease",
            "disease_status",
        ),
        "Ecotype": attributes.get("ecotype", ""),
        "collection date": attributes.get("collection_date", ""),
    }

rows = []

for _, run_info in runs.iterrows():
    biosample_id = run_info.get("BioSample", "").strip()
    run_id = run_info.get("Run", "").strip()
    sample = biosamples.get(biosample_id, {})

    rows.append({
        "Organism": run_info.get("ScientificName", ORGANISM),
        "Run Id": run_id,
        "Sample Accession": biosample_id,
        "Breed": sample.get("Breed", ""),
        "Links": f"https://www.ncbi.nlm.nih.gov/sra/?term={urlquote(run_id)}",
        "Experiment Accession": run_info.get("Experiment", ""),
        "BioProject": run_info.get("BioProject", ""),
        "Study": run_info.get("SRAStudy", ""),
        "Description": sample.get("Description", ""),
        "Geographical location": sample.get("Geographical location", ""),
        "Tissue": sample.get("Tissue", ""),
        "Sex": sample.get("Sex", ""),
        "Birth date": sample.get("Birth date", ""),
        "Death date": sample.get("Death date", ""),
        "Disease": sample.get("Disease", ""),
        "Ecotype": sample.get("Ecotype", ""),
        "collection date": sample.get("collection date", ""),
        "Library strategy": run_info.get("LibraryStrategy", ""),
        "Library selection": run_info.get("LibrarySelection", ""),
        "Library source": run_info.get("LibrarySource", ""),
        "Instrument model": run_info.get("Model", ""),
        "Platform": run_info.get("Platform", ""),
        "Bases": run_info.get("bases", ""),
        "Reads": run_info.get("spots", ""),
    })

columns = [
    "Organism", "Run Id", "Sample Accession", "Breed", "Links",
    "Experiment Accession", "BioProject", "Study", "Description",
    "Geographical location", "Tissue", "Sex", "Birth date", "Death date",
    "Disease", "Ecotype", "collection date", "Library strategy",
    "Library selection", "Library source", "Instrument model", "Platform",
    "Bases", "Bases (G/M/K)", "Reads", "Reads (G/M/K)",
]

df = pd.DataFrame(rows)

for column in ["Bases", "Reads"]:
    df[column] = pd.to_numeric(df[column], errors="coerce").astype("Int64")

df["Bases (G/M/K)"] = df["Bases"].apply(human_readable_count)
df["Reads (G/M/K)"] = df["Reads"].apply(human_readable_count)

df = df.reindex(columns=columns)

sheet_name = "Felis_silvestris_runs"
df.to_excel(OUTPUT_XLSX, index=False, sheet_name=sheet_name)

workbook = load_workbook(OUTPUT_XLSX)
worksheet = workbook[sheet_name]

worksheet.freeze_panes = "A2"
worksheet.auto_filter.ref = worksheet.dimensions

for cell in worksheet[1]:
    cell.font = Font(bold=True)

headers = [cell.value for cell in worksheet[1]]

hyperlink_columns = {
    "Run Id": lambda value: f"https://www.ncbi.nlm.nih.gov/sra/?term={urlquote(value)}",
    "Sample Accession": lambda value: f"https://www.ncbi.nlm.nih.gov/biosample/{urlquote(value)}",
    "Links": lambda value: value,
    "Experiment Accession": lambda value: f"https://www.ncbi.nlm.nih.gov/sra/?term={urlquote(value)}",
    "BioProject": lambda value: f"https://www.ncbi.nlm.nih.gov/bioproject/{urlquote(value)}",
    "Study": lambda value: f"https://www.ncbi.nlm.nih.gov/sra/?term={urlquote(value)}",
}

for column_name, make_url in hyperlink_columns.items():
    column_number = headers.index(column_name) + 1

    for row in range(2, worksheet.max_row + 1):
        cell = worksheet.cell(row=row, column=column_number)

        if cell.value:
            cell.hyperlink = make_url(str(cell.value))
            cell.font = Font(color="0563C1", underline="single")

for column_cells in worksheet.columns:
    width = min(
        max(len(str(cell.value or "")) for cell in column_cells) + 2,
        40,
    )
    worksheet.column_dimensions[
        get_column_letter(column_cells[0].column)
    ].width = width

for column_name in ["Bases", "Reads"]:
    column_number = headers.index(column_name) + 1

    for row in range(2, worksheet.max_row + 1):
        worksheet.cell(
            row=row,
            column=column_number,
        ).number_format = "#,##0"

workbook.save(OUTPUT_XLSX)

print(f"\nCreated: {OUTPUT_XLSX}")
print(f"Run metadata rows exported: {len(df)}")
print(f"Unique Run IDs: {df['Run Id'].nunique()}")
print(f"Unique BioSample accessions: {df['Sample Accession'].nunique()}")
PY
