#we are using Strelka and its germline pipeline 
#first we made a bam_path.txt file which has all the paths of all pf the bam files 
#command for installing strelka 
conda create -n strelka_env \
    --override-channels \
    -c conda-forge \
    -c bioconda \
    --strict-channel-priority \
    strelka=2.9.10 \
    python=2.7 \
    samtools \
    bcftools \
    htslib \
    -y

#Open a new script file
nano /home/aayushi/strelka_felis/run_strelka.sh

#entire script
#!/usr/bin/env bash
set -euo pipefail

BAM_LIST="/home/aayushi/strelka_felis/bam_paths.txt"
REF="/home/aayushi/data1/reference/GCF_027475565.1_VMU_Ajub_asm_v1.0_genomic.fasta"
OUT_DIR="/home/aayushi/strelka_felis/results"
JOBS=80
MEM_GB=140

mkdir -p "$OUT_DIR"

# Index the reference only if its index is missing.
if [[ ! -f "${REF}.fai" ]]; then
    samtools faidx "$REF"
fi

# Validate BAM paths and ensure output folder names are unique.
declare -A seen=()
bams=()

while IFS= read -r bam || [[ -n "$bam" ]]; do
    bam="${bam%$'\r'}"
    [[ -z "$bam" ]] && continue

    if [[ "$bam" != /* || ! -f "$bam" ]]; then
        echo "Invalid or missing absolute BAM path: $bam" >&2
        exit 1
    fi

    sample=$(basename "$bam" .bam)

    if [[ -n "${seen[$sample]:-}" ]]; then
        echo "Duplicate BAM basename: $sample" >&2
        exit 1
    fi

    seen["$sample"]=1
    bams+=("$bam")
done < "$BAM_LIST"

if (( ${#bams[@]} == 0 )); then
    echo "No BAM paths found." >&2
    exit 1
fi

echo "Found ${#bams[@]} BAM files"

# Process one BAM at a time.
for bam in "${bams[@]}"; do
    sample=$(basename "$bam" .bam)
    run_dir="$OUT_DIR/$sample"

    echo "Starting: $sample"

    # Index the BAM only if no supported index exists.
    if [[ ! -f "${bam}.bai" &&
          ! -f "${bam%.bam}.bai" &&
          ! -f "${bam}.csi" ]]; then
        samtools index "$bam"
    fi

    # Configure a new workflow, or reuse an existing one.
    if [[ ! -f "$run_dir/runWorkflow.py" ]]; then
        configureStrelkaGermlineWorkflow.py \
            --bam "$bam" \
            --referenceFasta "$REF" \
            --runDir "$run_dir"
    fi

    # Run with up to 80 parallel jobs and a 140 GB memory budget.
    "$run_dir/runWorkflow.py" \
        -m local \
        -j "$JOBS" \
        --memGb "$MEM_GB" \
        > "$run_dir/execution.log" 2>&1

    echo "Finished: $sample"
done

echo "All samples completed."

#Check the script’s syntax:
bash -n /home/aayushi/strelka_felis/run_strelka.sh

#start the script 
nohup bash /home/aayushi/strelka_felis/run_strelka.sh \
    > /home/aayushi/strelka_felis/strelka_batch.log 2>&1 &

#monitor the progress 
tail -f /home/aayushi/strelka_felis/strelka_batch.log

#then I faced an error while running a script 
cat /home/aayushi/strelka_felis/strelka_batch.log
#Invalid or missing absolute BAM path: /home/aayushi/merged/lybica/SAMEA112127404.GCF_027475565.1_VMU_Ajub_asm_v1.0.bam
#If the file is listed, back up your path list and remove trailing whitespace
cp /home/aayushi/strelka_felis/bam_paths.txt \
   /home/aayushi/strelka_felis/bam_paths.txt.bak

sed -i 's/[[:space:]]*$//' /home/aayushi/strelka_felis/bam_paths.txt
#Then restart:
nohup bash /home/aayushi/strelka_felis/run_strelka.sh \
    > /home/aayushi/strelka_felis/strelka_batch.log 2>&1 &

#To see how long your batch has been running
ps -p 1907213 -o pid,etime,cmd

#To inspect the current sample’s progress
tail -n 30 /home/aayushi/strelka_felis/results/SRR26919620.GCF_027475565.1_VMU_Ajub_asm_v1.0/execution.log
