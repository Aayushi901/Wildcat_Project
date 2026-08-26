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

# Another command for the same above
#!/bin/sh

THREADS=8
OUT="$(pwd)/bam_paths.txt"
MISSING="$(pwd)/missing_srr_drr_ids.txt"
IDS="$(mktemp)"
ALL_BAMS="$(mktemp)"

trap 'rm -f "$ALL_BAMS" "$IDS"' 0 1 2 15

: > "$OUT"
: > "$MISSING"

printf '%s\n' \
SRR6071639 \
DRR308102 \
SRR6071640 \
SRR28242800 \
SRR22528610 \
SRR22528611 \
SRR22528612 \
SRR22528613 \
SRR22528614 \
SRR22528615 \
SRR22528616 \
SRR22528617 \
SRR22528618 \
SRR22528619 \
SRR16086860 \
SRR12381567 \
SRR12381568 \
SRR12381569 \
SRR12381570 \
SRR6071638 \
SRR4426169 \
DRR487136 \
DRR487137 \
DRR487138 \
DRR404790 \
DRR404791 \
DRR404792 \
DRR404793 \
SRR27902933 \
SRR12914278 \
SRR12381822 \
SRR12381823 \
SRR12381824 \
SRR12381825 \
SRR4426180 \
SRR4426181 \
SRR4426182 \
SRR4426170 \
SRR4426171 \
SRR4426172 \
SRR4426173 \
SRR4426174 \
SRR4426175 \
SRR4426176 \
SRR4426177 \
SRR4426178 \
SRR4426162 \
SRR4426163 \
SRR4426164 \
SRR4426165 \
SRR4426166 \
SRR4426167 \
SRR3042235 \
SRR2062628 \
SRR28762828 \
SRR28762833 \
SRR28762843 \
SRR28762852 \
SRR28762864 \
SRR28762865 > "$IDS"

find / -maxdepth 1 -type f -iname '*.bam' -print 2>/dev/null > "$ALL_BAMS"

find / -mindepth 1 -maxdepth 1 -type d -print0 2>/dev/null |
xargs -0 -n 1 -P "$THREADS" sh -c \
'find "$1" -type f -iname "*.bam" -print 2>/dev/null' sh >> "$ALL_BAMS"

awk -F/ -v missing="$MISSING" '
NR == FNR {
  ids[++n] = $1
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
' "$IDS" "$ALL_BAMS" > "$OUT"

echo "Found BAM paths: $OUT"
echo "Missing IDs: $MISSING"
