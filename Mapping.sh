#command for having fastqc and multiqc 
conda create -n fastqc_env -c conda-forge -c bioconda fastqc multiqc
#command for doing fastqc 
find /home/aayushi/felis_catus/fastq_files -name "*.fastq.gz" | \
xargs -P 4 -n 1 fastqc -t 4 -o /home/aayushi/felis_catus/fastqc_results

