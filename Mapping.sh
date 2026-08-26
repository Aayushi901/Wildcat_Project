#command for having fastqc and multiqc 
conda create -n fastqc_env -c conda-forge -c bioconda fastqc multiqc
#command for doing fastqc 
find /home/aayushi/felis_catus/fastq_files -name "*.fastq.gz" | \
xargs -P 4 -n 1 fastqc -t 4 -o /home/aayushi/felis_catus/fastqc_results
#frm sra to fastgz
mkdir -p /home/aayushi/felis_catus/fastq_files/nap_sib_sia_toy
mkdir -p /home/aayushi/felis_catus/tmp

xargs -P 10 -n 1 bash -c '
echo "Processing $1"

fasterq-dump "$1" \
  --split-files \
  --threads 12 \
  --outdir /home/aayushi/felis_catus/fastq_files/nap_sib_sia_toy \
  --temp /home/aayushi/felis_catus/tmp &&

gzip /home/aayushi/felis_catus/fastq_files/nap_sib_sia_toy/${1}_*.fastq
' _ <<EOF
ERR3485781
ERR3493827
SRR11392454
SRR11392455
SRR11392605
SRR11392638
SRR13149608
SRR13149609
SRR13149610
SRR15512204
SRR15512212
SRR18908372
SRR5051110
SRR5363122
SRR5363127
SRR5363128
SRR5363129
SRR5366693
SRR5373737
SRR5373738
SRR5373739
EOF

#find the bam files throughout the server and list the missing bam files as well
OUT="$PWD/bam_paths.txt"
MISSING="$PWD/missing_srr_drr_ids.txt"

find / -type f -iname '*.bam' -print 2>/dev/null |
awk -F/ -v missing="$MISSING" '
  NR == FNR {
    ids[++n] = $1
    wanted[$1] = 1
    next
  }
  {
    filename = toupper($NF)
    for (i = 1; i <= n; i++) {
      id = ids[i]
      if (index(filename, id)) {
        print
        found[id] = 1
        break
      }
    }
  }
  END {
    for (i = 1; i <= n; i++)
      if (!found[ids[i]]) print ids[i] > missing
  }
' <(
  printf '%s\n' \
    SRR6071639 DRR308102 SRR6071640 SRR28242800 \
    SRR225286{10..19} SRR16086860 SRR123815{67..70} \
    SRR6071638 SRR4426169 DRR4871{36..38} DRR40479{0..3} \
    SRR27902933 SRR12914278 SRR123818{22..25} \
    SRR442618{0..2} SRR442617{0..8} SRR442616{2..7} \
    SRR3042235 SRR2062628 \
    SRR287628{28,33,43,52,64,65}
) - > "$OUT"

echo "Found BAM paths: $OUT"
echo "IDs without a BAM: $MISSING"
