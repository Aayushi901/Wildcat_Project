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
