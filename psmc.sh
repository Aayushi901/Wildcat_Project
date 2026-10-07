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
