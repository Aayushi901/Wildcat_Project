 conda create -n psmc-linux   --override-channels   -c conda-forge -c bioconda   --strict-channel-priority   samtools=1.24   bcftools=1.24   psmc=0.6.5   bedtools=2.31.1   perl   gnuplot
  979  conda activate psmc-linux
  980  mkdir -p /home/aayushi/psmc_linux
  981  cd /home/aayushi/psmc_linux
  982  conda list --explicit > software_versions.txt
  983  samtools --version | head -n 2
  984  bcftools --version | head -n 2
  985  command -v psmc
  986  command -v fq2psmcfa
  987  command -v splitfa
  988  command -v vcfutils.pl
  989  command -v psmc_plot.pl
  990  command -v gnuplot
  991  ls /home/aayushi/data/bam_files5/SRR6071639.GCF_027475565.1_VMU_Ajub_asm_v1.0.bam
  992  ls /home/aayushi/data/bam_files5/SRR6071640.GCF_027475565.1_VMU_Ajub_asm_v1.0.bam
  993  conda activate psmc-linux
  994  mkdir -p /home/aayushi/data1/psmc_linux
  995  cd /home/aayushi/data1/psmc_linux
  996  mkdir -p reference qc flat_headed rusty_spotted
  997  REF="/home/aayushi/data1/reference/GCF_027475565.1_VMU_Ajub_asm_v1.0_genomic.fasta"
  998  FLAT="/home/aayushi/data/bam_files5/SRR6071639.GCF_027475565.1_VMU_Ajub_asm_v1.0.bam"
  999  RUSTY="/home/aayushi/data/bam_files5/SRR6071640.GCF_027475565.1_VMU_Ajub_asm_v1.0.bam"
 1000  ls -lh "$REF" "$FLAT" "$RUSTY"
 1001  samtools quickcheck -v "$FLAT" "$RUSTY"   && echo "PASS: both BAM files passed quickcheck"
 1002  ls
 1003  samtools view -H "$FLAT" > qc/flat_headed.header.sam
 1004  samtools view -H "$RUSTY" > qc/rusty_spotted.header.sam
 1005  echo "=== FLAT-HEADED CAT ==="
 1006  awk '$1=="@HD" || $1=="@RG"' qc/flat_headed.header.sam
 1007  echo "=== RUSTY-SPOTTED CAT ==="
 1008  awk '$1=="@HD" || $1=="@RG"' qc/rusty_spotted.header.sam
 1009  samtools dict "$REF" > reference/cheetah.dict
 1010  for item in   reference/cheetah.dict   qc/flat_headed.header.sam   qc/rusty_spotted.header.sam; do   awk -F '\t' '
    BEGIN { OFS="\t" }
    $1=="@SQ" {
      name=""; len=""
      for (i=2; i<=NF; i++) {
        if ($i ~ /^SN:/) name=substr($i,4)
        if ($i ~ /^LN:/) len=substr($i,4)
      }
      print name,len
    }
  ' "$item" | LC_ALL=C sort -k1,1 > "${item}.names_lengths.tsv"; done
 1011  diff -u   reference/cheetah.dict.names_lengths.tsv   qc/flat_headed.header.sam.names_lengths.tsv   && echo "PASS: flat-headed BAM sequence names and lengths match"
 1012  diff -u   reference/cheetah.dict.names_lengths.tsv   qc/rusty_spotted.header.sam.names_lengths.tsv   && echo "PASS: rusty-spotted BAM sequence names and lengths match"
 1013  awk '
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
 1014  awk '
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
 1015  samtools flagstat -@ 4 "$FLAT" > qc/flat_headed.flagstat.txt
 1016  samtools flagstat -@ 4 "$RUSTY" > qc/rusty_spotted.flagstat.txt
 1017  echo "=== FLAT-HEADED CAT ==="
 1018  cat qc/flat_headed.flagstat.txt
 1019  echo "=== RUSTY-SPOTTED CAT ==="
 1020  cat qc/rusty_spotted.flagstat.txt
 1021  conda activate psmc-linux
 1022  cd /home/aayushi/data1/psmc_linux
 1023  conda list --explicit > software_versions.txt
 1024  echo "=== FLAT-HEADED CAT: PROCESSING HISTORY ==="
 1025  awk -F '\t' '$1=="@PG" || $1=="@CO"'   qc/flat_headed.header.sam
 1026  echo "=== RUSTY-SPOTTED CAT: PROCESSING HISTORY ==="
 1027  awk -F '\t' '$1=="@PG" || $1=="@CO"'   qc/rusty_spotted.header.sam
 1028  ls
 1029  cd flat_headed/
 1030  ls
 1031  cd
 1032  ls
 1033  cd data1
 1034  ls
 1035  cd psmc_linux/
 1036  ls
 1037  conda activate psmc-linux
 1038  cd /home/aayushi/data1/psmc_linux
 1039  set -o pipefail
 1040  FLAT="/home/aayushi/data/bam_files5/SRR6071639.GCF_027475565.1_VMU_Ajub_asm_v1.0.bam"
 1041  RUSTY="/home/aayushi/data/bam_files5/SRR6071640.GCF_027475565.1_VMU_Ajub_asm_v1.0.bam"
 1042  df -h .
 1043  cat > reference/autosomes.txt <<'EOF'
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

 1044  awk -F '\t' '
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
' reference/autosomes.txt reference/cheetah.dict   > reference/autosomes.bed
 1045  cat reference/autosomes.bed
 1046  awk '
{ total += $3-$2 }
END {
    print "Autosomes:", NR
    printf "Total autosomal bases: %.0f\n", total
}
' reference/autosomes.bed
 1047  samtools view -@ 4 -b   -q 30   -F 3844   -e '((flag & 1) == 0) || ((flag & 2) != 0)'   -L reference/autosomes.bed   -o flat_headed/autosomes.filtered.bam   "$FLAT"   && echo "PASS: flat-headed filtering completed"
 1048  samtools view -@ 4 -b   -q 30   -F 3844   -e '((flag & 1) == 0) || ((flag & 2) != 0)'   -L reference/autosomes.bed   -o rusty_spotted/autosomes.filtered.bam   "$RUSTY"   && echo "PASS: rusty-spotted filtering completed"
 1049  samtools quickcheck -v   flat_headed/autosomes.filtered.bam   rusty_spotted/autosomes.filtered.bam   && echo "PASS: both filtered BAMs passed quickcheck"
 1050  samtools index -@ 4 -c flat_headed/autosomes.filtered.bam   && echo "PASS: flat-headed BAM indexed"
 1051  samtools index -@ 4 -c rusty_spotted/autosomes.filtered.bam   && echo "PASS: rusty-spotted BAM indexed"
 1052  ( set -euo pipefail;  EXPECTED=$(awk '{n += $3-$2} END {printf "%.0f\n",n}' \
  reference/autosomes.bed);  for SAMPLE in flat_headed rusty_spotted; do     echo "Calculating coverage for: $SAMPLE";      samtools depth       -aa       -s       -q 20       -Q 30       -b reference/autosomes.bed       "$SAMPLE/autosomes.filtered.bam"     | awk       -v expected="$EXPECTED"       -v summary="qc/${SAMPLE}.depth_summary.txt"       -v hist="qc/${SAMPLE}.depth_histogram.tsv"       -v chrom="qc/${SAMPLE}.depth_by_chromosome.tsv" '
      BEGIN {
          OFS="\t"
          print "chromosome","positions","mean_depth" > chrom
      }
      {
          if ($1 != previous && previous != "") {
              print previous,nchr,sumchr/nchr > chrom
              nchr=0
              sumchr=0
          }

          previous=$1
          nchr++
          sumchr+=$3

          n++
          total+=$3
          histogram[$3]++

          if ($3==0) zero++
          if ($3>=10) ten++
          if ($3>=20) twenty++
      }
      END {
          if (n != expected) {
              printf "ERROR: expected %.0f positions, received %.0f\n",
                     expected,n > "/dev/stderr"
              exit 1
          }

          print previous,nchr,sumchr/nchr > chrom

          printf "Autosomal positions: %.0f\n",n > summary
          printf "Mean depth: %.3f\n",total/n > summary
          printf "Zero coverage: %.2f%%\n",100*zero/n > summary
          printf "Positions with depth >=10: %.2f%%\n",
                 100*ten/n > summary
          printf "Positions with depth >=20: %.2f%%\n",
                 100*twenty/n > summary

          for (d in histogram)
              printf "%d\t%.0f\n",d,histogram[d] > hist
      }
    ';      LC_ALL=C sort -k1,1n       "qc/${SAMPLE}.depth_histogram.tsv"       -o "qc/${SAMPLE}.depth_histogram.tsv";      echo "PASS: coverage completed for $SAMPLE"; done; )
 1053  echo "=== FLAT-HEADED CAT ==="
 1054  cat qc/flat_headed.depth_summary.txt
 1055  cat qc/flat_headed.depth_by_chromosome.tsv
 1056  echo "=== RUSTY-SPOTTED CAT ==="
 1057  cat qc/rusty_spotted.depth_summary.txt
 1058  cat qc/rusty_spotted.depth_by_chromosome.tsv
 1059  conda activate psmc-linux
 1060  cd /home/aayushi/data1/psmc_linux
 1061  awk -v low=10 -v high=36 '
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
 1062  awk -v low=10 -v high=27 '
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
 1063  REF="/home/aayushi/data1/reference/GCF_027475565.1_VMU_Ajub_asm_v1.0_genomic.fasta"
 1064  LC_ALL=C awk '
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
' reference/autosomes.txt "$REF"   > qc/reference_mask_check.txt   && cat qc/reference_mask_check.txt
 1065  ls /home/aayushi/data1/reference/
 1066  conda activate psmc-linux
 1067  cd /home/aayushi/data1/psmc_linux
 1068  set -o pipefail
 1069  REF="/home/aayushi/data1/reference/GCF_027475565.1_VMU_Ajub_asm_v1.0_genomic.fasta"
 1070  command -v bedtools
 1071  LC_ALL=C awk '
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
' reference/autosomes.txt "$REF"   | LC_ALL=C sort --parallel=8 -S 512M -T reference -k1,1 -k2,2n   | bedtools merge -i -   > reference/autosomes.softmask_and_ambiguity.bed.partial   && mv reference/autosomes.softmask_and_ambiguity.bed.partial         reference/autosomes.softmask_and_ambiguity.bed
 1072  ls /home/aayushi/data1/reference/
 1073  cd /home/aayushi/data1/psmc_linux
 1074  ls -lh reference/autosomes.softmask_and_ambiguity.bed
 1075  awk '
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
 1076  history
(psmc-linux) [aayushi@node02 psmc_linux]$ 
