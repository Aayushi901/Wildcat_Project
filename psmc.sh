#command for installing the essentials for doing PSMC
conda create -n psmc-linux \
  --override-channels \
  -c conda-forge -c bioconda \
  --strict-channel-priority \
  samtools=1.24 \
  bcftools=1.24 \
  psmc=0.6.5 \
  bedtools=2.31.1 \
  perl \
  gnuplot

#
conda activate psmc-linux

mkdir -p /home/aayushi/data1/psmc_linux
cd /home/aayushi/data1/psmc_linux

conda list --explicit > software_versions.txt

samtools --version | head -n 2
bcftools --version | head -n 2

command -v psmc
command -v fq2psmcfa
command -v splitfa
command -v vcfutils.pl
command -v psmc_plot.pl
command -v gnuplot

#providing the path of the files and creating directories 
  conda activate psmc-linux

mkdir -p /home/aayushi/data1/psmc_linux
cd /home/aayushi/data1/psmc_linux

mkdir -p reference qc flat-headed rusty-spotted

REF="/home/aayushi/data1/reference/GCF_027475565.1_VMU_Ajub_asm_v1.0_genomic.fasta"

FLAT="/home/aayushi/data/bam_files5/SRR6071639.GCF_027475565.1_VMU_Ajub_asm_v1.0.bam"

RUSTY="/home/aayushi/data/bam_files5/SRR6071640.GCF_027475565.1_VMU_Ajub_asm_v1.0.bam"

ls -lh "$REF" "$FLAT" "$RUSTY"

#Check whether the BAm headers and file endings are intact 
samtools quickcheck -v "$FLAT" "$RUSTY" \
  && echo "PASS: both BAM files passed quickcheck"

#Save the BAM headers and inspect sorting and sample information
  samtools view -H "$FLAT" > qc/flat_headed.header.sam

samtools view -H "$RUSTY" > qc/rusty_spotted.header.sam

echo "=== FLAT-HEADED CAT ==="
awk '$1=="@HD" || $1=="@RG"' qc/flat_headed.header.sam

echo "=== RUSTY-SPOTTED CAT ==="
awk '$1=="@HD" || $1=="@RG"' qc/rusty_spotted.header.sam

#Build a reference dictionary
#This reads your FASTA and writes a small dictionary inside the new working directory.
samtools dict "$REF" > reference/cheetah.dict


#Extract sequence names and lengths from the reference and both BAM headers:
for item in \
  reference/cheetah.dict \
  qc/flat_headed.header.sam \
  qc/rusty_spotted.header.sam
do
  awk -F '\t' '
    BEGIN { OFS="\t" }
    $1=="@SQ" {
      name=""; len=""
      for (i=2; i<=NF; i++) {
        if ($i ~ /^SN:/) name=substr($i,4)
        if ($i ~ /^LN:/) len=substr($i,4)
      }
      print name,len
    }
  ' "$item" | LC_ALL=C sort -k1,1 > "${item}.names_lengths.tsv"
done


#Compare each BAM with the reference:
diff -u \
  reference/cheetah.dict.names_lengths.tsv \
  qc/flat_headed.header.sam.names_lengths.tsv \
  && echo "PASS: flat-headed BAM sequence names and lengths match"

  diff -u \
  reference/cheetah.dict.names_lengths.tsv \
  qc/rusty_spotted.header.sam.names_lengths.tsv \
  && echo "PASS: rusty-spotted BAM sequence names and lengths match"

#Matching names and lengths establishes compatibility, but does not alone prove identical reference bases. Check whether the BAM headers contain reference checksums:
awk '
  $1=="@SQ" {
    total++
    for (i=2; i<=NF; i++)
      if ($i ~ /^M5:/) checksums++
  }
  END {
    print FILENAME ": reference sequences =", total,
          "; sequences with MD5 checksums =", checksums+0
  }
' qc/flat_headed.header.sam

awk '
  $1=="@SQ" {
    total++
    for (i=2; i<=NF; i++)
      if ($i ~ /^M5:/) checksums++
  }
  END {
    print FILENAME ": reference sequences =", total,
          "; sequences with MD5 checksums =", checksums+0
  }
' qc/rusty_spotted.header.sam

#Obtain mapping, pairing, and duplicate statistics
#These statistics report mapping, paired-read, and duplicate flags; they do not measure sequencing depth.
samtools flagstat -@ 4 "$FLAT" > qc/flat_headed.flagstat.txt

samtools flagstat -@ 4 "$RUSTY" > qc/rusty_spotted.flagstat.txt


echo "=== FLAT-HEADED CAT ==="
cat qc/flat_headed.flagstat.txt

echo "=== RUSTY-SPOTTED CAT ==="
cat qc/rusty_spotted.flagstat.txt

#
conda activate psmc-linux
cd /home/aayushi/data1/psmc_linux

conda list --explicit > software_versions.txt

echo "=== FLAT-HEADED CAT: PROCESSING HISTORY ==="
awk -F '\t' '$1=="@PG" || $1=="@CO"' \
  qc/flat_headed.header.sam

echo "=== RUSTY-SPOTTED CAT: PROCESSING HISTORY ==="
awk -F '\t' '$1=="@PG" || $1=="@CO"' \
  qc/rusty_spotted.header.sam
  
#For the next step, we will retain:
#High-quality reads flagged as unpaired, which include collapsed reads and singletons.
#High-quality paired reads only when they are properly paired.
#1. Set the paths
conda activate psmc-linux
cd /home/aayushi/data1/psmc_linux

set -o pipefail

FLAT="/home/aayushi/data/bam_files5/SRR6071639.GCF_027475565.1_VMU_Ajub_asm_v1.0.bam"

RUSTY="/home/aayushi/data/bam_files5/SRR6071640.GCF_027475565.1_VMU_Ajub_asm_v1.0.bam"

df -h .
#2. Select the 18 autosomes
#These accessions come from the cheetah assembly report. Chromosome X (NC_069389.1) is excluded for both cats.
cat > reference/autosomes.txt <<'EOF'
NC_069380.1
NC_069381.1
NC_069382.1
NC_069383.1
NC_069384.1
NC_069385.1
NC_069386.1
NC_069387.1
NC_069388.1
NC_069390.1
NC_069391.1
NC_069392.1
NC_069393.1
NC_069394.1
NC_069395.1
NC_069396.1
NC_069397.1
NC_069398.1
EOF

#Generate full-chromosome intervals using the lengths
awk -F '\t' '
BEGIN { OFS="\t" }
NR==FNR { wanted[$1]=1; next }
$1=="@SQ" {
    name=""; len=""
    for (i=2; i<=NF; i++) {
        if ($i ~ /^SN:/) name=substr($i,4)
        if ($i ~ /^LN:/) len=substr($i,4)
    }
    if (name in wanted) {
        print name,0,len
        found++
    }
}
END {
    if (found!=18) {
        print "ERROR: expected 18 autosomes, found " found > "/dev/stderr"
        exit 1
    }
}
' reference/autosomes.txt reference/cheetah.dict \
  > reference/autosomes.bed
#Check the result:
cat reference/autosomes.bed

awk '
{ total += $3-$2 }
END {
    print "Autosomes:", NR
    printf "Total autosomal bases: %.0f\n", total
}
' reference/autosomes.bed

#The filters require mapping quality ≥30, exclude unmapped/secondary/supplementary/QC-failed/duplicate-flagged records, and retain unpaired reads or properly paired reads. They select autosomal alignments only
samtools view -@ 4 -b \
  -q 30 \
  -F 3844 \
  -e '((flag & 1) == 0) || ((flag & 2) != 0)' \
  -L reference/autosomes.bed \
  -o flat_headed/autosomes.filtered.bam \
  "$FLAT" \
  && echo "PASS: flat-headed filtering completed"
  #Then the rusty-spotted cat:
  samtools view -@ 4 -b \
  -q 30 \
  -F 3844 \
  -e '((flag & 1) == 0) || ((flag & 2) != 0)' \
  -L reference/autosomes.bed \
  -o rusty_spotted/autosomes.filtered.bam \
  "$RUSTY" \
  && echo "PASS: rusty-spotted filtering completed"

  #Filtering preserves the original record order; the indexing step below will also detect incompatible sorting.
  samtools quickcheck -v \
  flat_headed/autosomes.filtered.bam \
  rusty_spotted/autosomes.filtered.bam \
  && echo "PASS: filtered BAMs passed quickcheck"
#
samtools index -@ 4 -c flat_headed/autosomes.filtered.bam \
  && echo "PASS: flat-headed BAM indexed"
#
samtools index -@ 4 -c rusty_spotted/autosomes.filtered.bam \
  && echo "PASS: rusty-spotted BAM indexed"

#1. Check how much sequence would remain under candidate depth filters
conda activate psmc-linux
cd /home/aayushi/data1/psmc_linux

awk -v low=10 -v high=36 '
{
    total += $2
    if ($1 < low) below += $2
    else if ($1 > high) above += $2
    else retained += $2
}
END {
    if (total==0) {
        print "ERROR: empty depth histogram"
        exit 1
    }
    print "FLAT-HEADED CAT"
    printf "Candidate depth range: %d-%d\n",low,high
    printf "Below range: %.2f%%\n",100*below/total
    printf "Above range: %.2f%%\n",100*above/total
    printf "Within range: %.2f%%\n",100*retained/total
}
' qc/flat_headed.depth_histogram.tsv

#
awk -v low=10 -v high=27 '
{
    total += $2
    if ($1 < low) below += $2
    else if ($1 > high) above += $2
    else retained += $2
}
END {
    if (total==0) {
        print "ERROR: empty depth histogram"
        exit 1
    }
    print "RUSTY-SPOTTED CAT"
    printf "Candidate depth range: %d-%d\n",low,high
    printf "Below range: %.2f%%\n",100*below/total
    printf "Above range: %.2f%%\n",100*above/total
    printf "Within range: %.2f%%\n",100*retained/total
}
' qc/rusty_spotted.depth_histogram.tsv

#2. Check whether your reference contains lowercase masking
REF="/home/aayushi/data1/reference/GCF_027475565.1_VMU_Ajub_asm_v1.0_genomic.fasta"

LC_ALL=C awk '
NR==FNR {
    wanted[$1]=1
    next
}
/^>/ {
    name=$1
    sub(/^>/,"",name)
    selected=(name in wanted)
    if (selected) chromosomes++
    next
}
selected {
    seq=$0
    gsub(/[[:space:]]/,"",seq)
    total += length(seq)

    copy=seq
    lowercase += gsub(/[acgt]/,"",copy)

    copy=seq
    ambiguous += gsub(/[^ACGTacgt]/,"",copy)
}
END {
    if (chromosomes!=18 || total!=2246314264) {
        print "ERROR: unexpected autosome count or sequence length" > "/dev/stderr"
        print "Chromosomes:",chromosomes > "/dev/stderr"
        printf "Bases: %.0f\n",total > "/dev/stderr"
        exit 1
    }

    print "Autosomes:",chromosomes
    printf "Total bases: %.0f\n",total
    printf "Lowercase A/C/G/T: %.0f (%.2f%%)\n",
           lowercase,100*lowercase/total
    printf "Ambiguous bases, including N: %.0f (%.2f%%)\n",
           ambiguous,100*ambiguous/total
}
' reference/autosomes.txt "$REF" \
  > qc/reference_mask_check.txt \
  && cat qc/reference_mask_check.txt

#The checks passed. reference contains 31.99% lowercase sequence, and the candidate depth filters retain 54.89% of the flat-headed genome and 67.19% of the rusty-spotted genome.
#These percentages overlap: we cannot subtract 31.99% from the depth-retained percentages to calculate usable sequence. Some lowercase regions already fail the depth filters.
#NCBI distributes genomes with lowercase masking from WindowMasker or RepeatMasker. Your result is consistent with that convention, but the lowercase count alone cannot identify the masking method.

#Next, extract the lowercase and ambiguous regions into a candidate exclusion BED.
#1. Set the paths
conda activate psmc-linux
cd /home/aayushi/data1/psmc_linux

set -o pipefail

REF="/home/aayushi/data1/reference/GCF_027475565.1_VMU_Ajub_asm_v1.0_genomic.fasta"

command -v bedtools
#2. Extract and merge the masked intervals
LC_ALL=C awk '
BEGIN { OFS="\t" }

NR==FNR {
    wanted[$1]=1
    next
}

/^>/ {
    name=$1
    sub(/^>/,"",name)
    selected=(name in wanted)
    position=0
    if (selected) chromosomes++
    next
}

selected {
    seq=$0
    gsub(/[[:space:]]/,"",seq)

    line_length=length(seq)
    total+=line_length
    offset=0

    while (match(seq,/[^ACGT]+/)) {
        start=position+offset+RSTART-1
        end=start+RLENGTH
        print name,start,end

        consumed=RSTART+RLENGTH-1
        offset+=consumed
        seq=substr(seq,consumed+1)
    }

    position+=line_length
}

END {
    if (chromosomes!=18 || total!=2246314264) {
        print "ERROR: unexpected autosome count or sequence length" > "/dev/stderr"
        exit 1
    }
}
' reference/autosomes.txt "$REF" \
  | LC_ALL=C sort --parallel=8 -S 512M -T reference -k1,1 -k2,2n \
  | bedtools merge -i - \
  > reference/autosomes.softmask_and_ambiguity.bed.partial \
  && mv reference/autosomes.softmask_and_ambiguity.bed.partial \
        reference/autosomes.softmask_and_ambiguity.bed

#This produces zero-based, end-exclusive BED coordinates. Merging joins adjacent intervals, including regions split across FASTA lines.
#Your previous counts imply:
#- Lowercase bases: 718,589,567
#- Ambiguous bases: 2,800
#- Total excluded bases: 718,592,367
#This command checks that the BED represents exactly that total:
awk '
{
    if ($2<0 || $3<=$2) bad++
    masked += $3-$2
}
END {
    printf "Mask intervals: %.0f\n",NR
    printf "Masked bases: %.0f\n",masked
    printf "Masked percentage: %.4f%%\n",
           100*masked/2246314264
    printf "Bases outside this mask: %.0f\n",
           2246314264-masked

    if (bad || masked!=718592367) {
        print "ERROR: mask validation failed"
        exit 1
    }

    print "PASS: mask length matches the FASTA counts"
}
' reference/autosomes.softmask_and_ambiguity.bed

#(
set -euo pipefail

REFDIR="/home/aayushi/data1/reference"
PREFIX="GCF_027475565.1_VMU_Ajub_asm_v1.0"

echo "Checking compressed file integrity..."
gzip -t "$REFDIR/${PREFIX}_genomic.fna.gz"

echo "Checksum of decompressed .fna.gz:"
gzip -dc "$REFDIR/${PREFIX}_genomic.fna.gz" | sha256sum

echo "Checksum of .fasta:"
sha256sum "$REFDIR/${PREFIX}_genomic.fasta"
)
#it shows
#Mask intervals: 12765917
#Masked bases: 718592367
#Masked percentage: 31.9898%
#Bases outside this mask: 1527721897
#PASS: mask length matches the FASTA counts

#next step- This step creates intermediate BCF files containing both variant and nonvariant calls at covered positions. Depth filtering and the exclusion mask will be applied during consensus preparation. The mask’s extraction is validated; its download/masking provenance still needs confirmation.
conda activate psmc-linux
cd /home/aayushi/data1/psmc_linux

(
set -euo pipefail

REF="/home/aayushi/data1/reference/GCF_027475565.1_VMU_Ajub_asm_v1.0_genomic.fasta"

mkdir -p flat_headed/bcf rusty_spotted/bcf logs/calling

if [ ! -e reference/cheetah.fa ]; then
    ln -s "$REF" reference/cheetah.fa
fi

samtools faidx reference/cheetah.fa

test -s flat_headed/autosomes.filtered.bam.csi
test -s rusty_spotted/autosomes.filtered.bam.csi

awk '{
    print "flat_headed", $1
    print "rusty_spotted", $1
}' reference/autosomes.txt > qc/calling_jobs.txt

echo "Number of chromosome jobs:"
wc -l < qc/calling_jobs.txt
)

#
xargs -P 8 -n 2 bash -c '
set -euo pipefail

SAMPLE="$1"
CHR="$2"

OUT="$SAMPLE/bcf/$CHR.bcf"
TMP="$SAMPLE/bcf/$CHR.partial.bcf"
LOG="logs/calling/${SAMPLE}.${CHR}"

if [ -s "$OUT" ] && [ -s "$OUT.csi" ] && [ -f "$OUT.done" ]; then
    echo "SKIP completed: $SAMPLE $CHR"
    exit 0
fi

echo "START: $SAMPLE $CHR"

bcftools mpileup \
    --full-BAQ \
    -E \
    -q 30 \
    -Q 20 \
    -d 10000 \
    -f reference/cheetah.fa \
    -r "$CHR" \
    -Ou \
    "$SAMPLE/autosomes.filtered.bam" \
    2> "${LOG}.mpileup.log" \
  | bcftools call \
      -c \
      --ploidy 2 \
      -Ob \
      -o "$TMP" \
      2> "${LOG}.call.log"

bcftools index -f "$TMP" 2> "${LOG}.index.log"

mv "$TMP" "$OUT"
mv "$TMP.csi" "$OUT.csi"
touch "$OUT.done"

echo "DONE: $SAMPLE $CHR"
' _ < qc/calling_jobs.txt \
  && echo "PASS: all 36 chromosome-calling jobs completed"

#
echo "Flat-headed completed chromosomes:"
find flat_headed/bcf -maxdepth 1 -name '*.bcf.done' | wc -l

echo "Rusty-spotted completed chromosomes:"
find rusty_spotted/bcf -maxdepth 1 -name '*.bcf.done' | wc -l

#Then check the sample identities:
bcftools query -l flat_headed/bcf/NC_069380.1.bcf

bcftools query -l rusty_spotted/bcf/NC_069380.1.bcf

#Next, convert the BCF files into intermediate diploid FASTQ files, eight chromosomes at a time. We will use the candidate depth ranges:
#1. Prepare the output folders
conda activate psmc-linux
cd /home/aayushi/data1/psmc_linux

mkdir -p \
  flat_headed/consensus_d10_D36 \
  rusty_spotted/consensus_d10_D27 \
  logs/consensus

#2. Run the conversion
xargs -P 8 -n 2 bash -c '
set -euo pipefail

SAMPLE="$1"
CHR="$2"

case "$SAMPLE" in
    flat_headed)
        MAXDP=36
        EXPECTED="SRR6071639"
        ;;
    rusty_spotted)
        MAXDP=27
        EXPECTED="SRR6071640"
        ;;
    *)
        echo "ERROR: unexpected sample $SAMPLE" >&2
        exit 1
        ;;
esac

INPUT="$SAMPLE/bcf/$CHR.bcf"
DIR="$SAMPLE/consensus_d10_D${MAXDP}"
OUT="$DIR/$CHR.raw.fq.gz"
LOG="logs/consensus/${SAMPLE}.${CHR}.d10_D${MAXDP}"

if [ -s "$OUT" ] && [ -f "$OUT.done" ]; then
    echo "SKIP completed: $SAMPLE $CHR"
    exit 0
fi

test -s "$INPUT"
test -f "$INPUT.done"

ACTUAL=$(bcftools query -l "$INPUT" 2> "${LOG}.sample.log")

if [ "$ACTUAL" != "$EXPECTED" ]; then
    echo "ERROR: unexpected sample identity in $INPUT" >&2
    exit 1
fi

bcftools view -h "$INPUT" \
    > "${LOG}.header.txt" \
    2> "${LOG}.header.log"

awk "
    /^##INFO=<ID=FQ,/ { fq=1 }
    /^##INFO=<ID=DP,/ { dp=1 }
    /^##INFO=<ID=MQ,/ { mq=1 }
    END { exit !(fq && dp && mq) }
" "${LOG}.header.txt" || {
    echo "ERROR: required FQ, DP or MQ annotation missing: $INPUT" >&2
    exit 1
}

echo "START: $SAMPLE $CHR"

bcftools view -Ov "$INPUT" \
    2> "${LOG}.view.log" \
  | vcfutils.pl vcf2fq \
      -d 10 \
      -D "$MAXDP" \
      -Q 30 \
      -l 10 \
      2> "${LOG}.vcf2fq.log" \
  | gzip -c > "$OUT.partial"

gzip -t "$OUT.partial"

mv "$OUT.partial" "$OUT"
touch "$OUT.done"

echo "DONE: $SAMPLE $CHR"
' _ < qc/calling_jobs.txt \
  && echo "PASS: all 36 intermediate FASTQ conversions completed"

#Here, vcf2fq -Q 30 filters the site-level mapping-quality annotation, and -l 10 marks the neighbourhood of indels. The converter represents filtered positions using lowercase sequence or missing bases; we will handle those explicitly in the next step.
echo "Flat-headed FASTQ files completed:"
find flat_headed/consensus_d10_D36 \
  -maxdepth 1 -name '*.raw.fq.gz.done' | wc -l

echo "Rusty-spotted FASTQ files completed:"
find rusty_spotted/consensus_d10_D27 \
  -maxdepth 1 -name '*.raw.fq.gz.done' | wc -l

  ##Both samples are ready for masking. We will now create chromosome FASTAs that:
#Retain consensus bases with quality ≥20.
#Replace lowercase calls and reference-mask regions with N.
#Preserve heterozygous bases as IUPAC codes.
#Fill any missing chromosome ends with N, preserving reference coordinates.
conda activate psmc-linux
cd /home/aayushi/data1/psmc_linux

conda install -y --freeze-installed \
  --override-channels -c conda-forge -c bioconda seqtk

conda list --explicit > software_versions.txt

#2. Prepare chromosome masks and job instructions
mkdir -p \
  reference/masks_by_chr \
  flat_headed/masked_d10_D36 \
  rusty_spotted/masked_d10_D27 \
  logs/masking

  #Split the existing mask into chromosome files:
  awk '
{
    file="reference/masks_by_chr/" $1 ".bed"
    if (file != previous) {
        if (previous != "") close(previous)
        printf "%s","" > file
        previous=file
    }
    print $0 > file
}
' reference/autosomes.softmask_and_ambiguity.bed

awk '
BEGIN { OFS="\t" }
NR==FNR { sizes[$1]=$3; next }
{
    if (!($2 in sizes)) {
        print "ERROR: missing chromosome length" > "/dev/stderr"
        exit 1
    }
    print $1,$2,sizes[$2]
}
' reference/autosomes.bed qc/calling_jobs.txt \
  > qc/masking_jobs.txt

  #3. Create the length-checking and reporting instructions
  #This is a small AWK helper used by the Linux commands, not Python. Copy the entire block, including the final AWK line:
  cat > qc/check_masked_fasta.awk <<'AWK'
BEGIN {
    OFS="\t"
    nblock=sprintf("%80s","")
    gsub(/ /,"N",nblock)
}

/^>/ {
    headers++
    name=$1
    sub(/^>/,"",name)

    if (headers!=1 || name!=chr) {
        print "ERROR: unexpected FASTA header" > "/dev/stderr"
        failed=1
        exit 1
    }

    print ">" chr
    next
}

{
    if (headers!=1) {
        print "ERROR: sequence without header" > "/dev/stderr"
        failed=1
        exit 1
    }

    seq=$0
    gsub(/[[:space:]]/,"",seq)

    # Keep only uppercase homozygous or biallelic heterozygous calls.
    # Lowercase calls, Ns, and other ambiguity codes become N.
    gsub(/[^ACGTMRWSYK]/,"N",seq)

    observed+=length(seq)
    if (observed>expected) {
        print "ERROR: sequence exceeds reference length" > "/dev/stderr"
        failed=1
        exit 1
    }

    copy=seq
    callable+=gsub(/[ACGTMRWSYK]/,"",copy)

    copy=seq
    heterozygous+=gsub(/[MRWSYK]/,"",copy)

    print seq
}

END {
    if (failed) exit 1

    if (headers!=1 || observed==0) {
        print "ERROR: empty or invalid FASTA" > "/dev/stderr"
        exit 1
    }

    # Complete any uncovered chromosome tail with missing bases.
    remaining=expected-observed
    while (remaining>0) {
        chunk=(remaining>80 ? 80 : remaining)
        print substr(nblock,1,chunk)
        remaining-=chunk
    }

    printf "%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\n",
           chr,expected,callable,heterozygous,
           expected-callable,expected-observed > report
}
AWK

#4. Apply the masking in parallel
xargs -P 8 -n 3 bash -c '
set -euo pipefail

SAMPLE="$1"
CHR="$2"
LENGTH="$3"

case "$SAMPLE" in
    flat_headed) MAXDP=36 ;;
    rusty_spotted) MAXDP=27 ;;
    *) echo "ERROR: unexpected sample" >&2; exit 1 ;;
esac

INPUT="$SAMPLE/consensus_d10_D${MAXDP}/$CHR.raw.fq.gz"
MASK="reference/masks_by_chr/$CHR.bed"
DIR="$SAMPLE/masked_d10_D${MAXDP}"
OUT="$DIR/$CHR.masked.fa.gz"
REPORT="$DIR/$CHR.qc.tsv"
LOG="logs/masking/${SAMPLE}.${CHR}"

if [ -s "$OUT" ] && [ -s "$REPORT" ] && [ -f "$OUT.done" ]; then
    echo "SKIP completed: $SAMPLE $CHR"
    exit 0
fi

test -s "$INPUT"
test -f "$INPUT.done"
test -s "$MASK"

echo "START: $SAMPLE $CHR"

seqtk seq \
    -q 20 \
    -Q 33 \
    -n N \
    -M "$MASK" \
    -A \
    -l 80 \
    "$INPUT" \
    2> "${LOG}.seqtk.log" \
  | LC_ALL=C awk \
      -v chr="$CHR" \
      -v expected="$LENGTH" \
      -v report="$REPORT.partial" \
      -f qc/check_masked_fasta.awk \
      2> "${LOG}.validation.log" \
  | gzip -c > "$OUT.partial"

gzip -t "$OUT.partial"
test -s "$REPORT.partial"

mv "$OUT.partial" "$OUT"
mv "$REPORT.partial" "$REPORT"
touch "$OUT.done"

echo "DONE: $SAMPLE $CHR"
' _ < qc/masking_jobs.txt \
  && echo "PASS: all 36 masked chromosome FASTAs completed"

#5. Display the usable-sequence summaries
for SAMPLE in flat_headed rusty_spotted
do
    case "$SAMPLE" in
        flat_headed) DIR="$SAMPLE/masked_d10_D36" ;;
        rusty_spotted) DIR="$SAMPLE/masked_d10_D27" ;;
    esac

    echo "=== $SAMPLE ==="

    echo "Completed chromosomes:"
    find "$DIR" -maxdepth 1 -name '*.masked.fa.gz.done' | wc -l

    awk '
    {
        chromosomes++
        total+=$2
        callable+=$3
        het+=$4
        missing+=$5
        tail+=$6
    }
    END {
        if (chromosomes!=18 || total!=2246314264) {
            print "ERROR: incomplete chromosome reports"
            exit 1
        }

        printf "Reference bases: %.0f\n",total
        printf "Retained callable bases: %.0f (%.2f%% of autosomes)\n",
               callable,100*callable/total
        printf "Missing or masked bases: %.0f (%.2f%%)\n",
               missing,100*missing/total
        printf "Callable fraction outside reference mask: %.2f%%\n",
               100*callable/1527721897
        printf "Retained heterozygous bases: %.0f\n",het
        if (callable>0)
            printf "Heterozygous fraction of callable bases: %.6f\n",
                   het/callable
        printf "Missing chromosome-tail bases added: %.0f\n",tail
    }
    ' "$DIR"/*.qc.tsv
done

#Both masked genomes are complete. We retained approximately 851 million callable bases for the flat-headed cat and 1.091 billion for the rusty-spotted cat.
#The missing fraction outside the reference mask is still 44.31% and 28.61%, respectively. We can proceed with exploratory PSMC runs, but the unequal coverage and filtering mean these heterozygosity estimates should not yet be interpreted as a biological difference between species.
#Next, create and validate one PSMC input per cat. Both conversions will run simultaneously.
#PSMC uses 100-base bins. Its converter can omit chromosomes with fewer than 20% usable bases; we will explicitly check that all 18 survive. We will also pad the final incomplete bin with N so it is treated consistently.
#1. Prepare the folders
conda activate psmc-linux
cd /home/aayushi/data1/psmc_linux

mkdir -p \
  flat_headed/psmc_d10_D36 \
  rusty_spotted/psmc_d10_D27 \
  logs/psmc_input

#2. Create the chromosome-length check and padding instructions
cat > qc/pad_psmc_input.awk <<'AWK'
function finish_chromosome(    padding) {
    if (chr=="") return

    if (observed!=sizes[chr]) {
        print "ERROR: incorrect sequence length for " chr > "/dev/stderr"
        failed=1
        exit 1
    }

    padding=(100-observed%100)%100
    if (padding>0) print substr(nblock,1,padding)
}

BEGIN {
    nblock=sprintf("%100s","")
    gsub(/ /,"N",nblock)
}

NR==FNR {
    sizes[$1]=$3
    next
}

/^>/ {
    finish_chromosome()

    chr=$1
    sub(/^>/,"",chr)

    if (!(chr in sizes) || seen[chr]++) {
        print "ERROR: unexpected or repeated chromosome " chr > "/dev/stderr"
        failed=1
        exit 1
    }

    chromosomes++
    observed=0
    print ">" chr
    next
}

{
    seq=$0
    gsub(/[[:space:]]/,"",seq)

    if (chr=="" || seq ~ /[^ACGTMRWSYKN]/) {
        print "ERROR: invalid masked FASTA" > "/dev/stderr"
        failed=1
        exit 1
    }

    observed+=length(seq)
    print seq
}

END {
    if (failed) exit 1
    finish_chromosome()

    if (chromosomes!=18) {
        print "ERROR: expected 18 chromosomes" > "/dev/stderr"
        exit 1
    }
}
AWK

#3. Convert both samples in parallel
printf '%s\n' flat_headed rusty_spotted \
  | xargs -P 2 -n 1 bash -c '
set -euo pipefail

SAMPLE="$1"

case "$SAMPLE" in
    flat_headed) MAXDP=36 ;;
    rusty_spotted) MAXDP=27 ;;
    *) exit 1 ;;
esac

INDIR="$SAMPLE/masked_d10_D${MAXDP}"
OUTDIR="$SAMPLE/psmc_d10_D${MAXDP}"
OUT="$OUTDIR/input.psmcfa"

echo "START: $SAMPLE"

gzip -dc "$INDIR"/*.masked.fa.gz \
  | LC_ALL=C awk \
      -f qc/pad_psmc_input.awk \
      reference/autosomes.bed - \
  | fq2psmcfa -s 100 - \
      > "$OUT.partial" \
      2> "logs/psmc_input/${SAMPLE}.log"

test -s "$OUT.partial"
mv "$OUT.partial" "$OUT"

echo "DONE: $SAMPLE"
' _ \
  && echo "PASS: both PSMC input conversions completed"

  #4. Validate both PSMC inputs
  The symbols mean:
#K: a retained bin containing a heterozygous call.
#T: a retained bin with no heterozygous call.
##N: a missing bin.
#Bin-level missingness differs from base-level missingness: the standard converter labels a bin missing when more than 90% of its bases are missing. Consequently, retained bins can still contain masked bases.
cat > qc/check_psmc_input.awk <<'AWK'
function check_length() {
    if (chr!="" && bins!=int((sizes[chr]+99)/100)) {
        print "ERROR: incorrect bin count for " chr > "/dev/stderr"
        bad=1
    }
}

NR==FNR {
    sizes[$1]=$3
    next
}

/^>/ {
    check_length()

    chr=$1
    sub(/^>/,"",chr)
    bins=0

    if (!(chr in sizes) || seen[chr]++) bad=1
    chromosomes++
    next
}

{
    seq=$0
    gsub(/[[:space:]]/,"",seq)

    if (chr=="" || seq ~ /[^KNT]/) bad=1

    bins+=length(seq)
    total+=length(seq)

    copy=seq
    het+=gsub(/K/,"",copy)
    copy=seq
    hom+=gsub(/T/,"",copy)
    copy=seq
    missing+=gsub(/N/,"",copy)
}

END {
    check_length()

    print "Chromosomes:",chromosomes
    printf "Total bins: %.0f\n",total
    printf "Heterozygous bins (K): %.0f\n",het
    printf "Nonheterozygous bins (T): %.0f\n",hom
    printf "Missing bins (N): %.0f\n",missing

    if (total>0)
        printf "Missing-bin percentage: %.2f%%\n",100*missing/total

    if (bad || chromosomes!=18 || het==0 || hom==0) {
        print "ERROR: PSMC input needs review"
        exit 1
    }

    print "PASS: chromosome identities, bin lengths and symbols validated"
}
AWK

#
echo "=== FLAT-HEADED CAT ==="

awk -f qc/check_psmc_input.awk \
  reference/autosomes.bed \
  flat_headed/psmc_d10_D36/input.psmcfa

#
echo "=== RUSTY-SPOTTED CAT ==="

awk -f qc/check_psmc_input.awk \
  reference/autosomes.bed \
  rusty_spotted/psmc_d10_D27/input.psmcfa
