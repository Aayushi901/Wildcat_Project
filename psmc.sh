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


  
