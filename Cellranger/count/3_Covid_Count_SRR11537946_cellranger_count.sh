#! /usr/bin/env bash

#BSUB -J CR3.0.0
#BSUB -n 56
#BSUB -R "span[hosts=1]"
#BSUB -W 999:59
#BSUB -o Covid_Count_SRR11537946_Cov3%J.out
#BSUB -e Covid_Count_SRR11537946_Cov3%J.err
#BSUB -N


/home/ajaffer/cellranger_3.0.1/cellranger-3.0.1/cellranger count --id=SRR11537946_Cov3 \
                 --transcriptome=/home/ajaffer/COVID_GEO_GSE145926/Cov_genome_3/refdata-cellranger-GRCh38-3.0.0_COVID \
                 --fastqs=/home/ajaffer/COVID_GEO_GSE145926/FASTQ_2 \
                 --sample=SRR11537946 \
                 --localcores=56

