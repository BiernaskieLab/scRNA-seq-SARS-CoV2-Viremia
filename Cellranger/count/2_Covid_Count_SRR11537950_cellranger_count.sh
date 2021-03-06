#! /usr/bin/env bash

#BSUB -J CR3.0.0
#BSUB -n 56
#BSUB -R "span[hosts=1]"
#BSUB -W 999:59
#BSUB -o Covid_Count_SRR11537950_Cov%J.out
#BSUB -e Covid_Count_SRR11537950_Cov%J.err
#BSUB -N


/home/ajaffer/cellranger_3.0.1/cellranger-3.0.1/cellranger count --id=SRR11537950_Cov \
                 --transcriptome=/home/ajaffer/COVID_GEO_GSE145926/Cov_genome_2/mkref/refdata-cellranger-GRCh38-3.0.0_COVID \
                 --fastqs=/home/ajaffer/COVID_GEO_GSE145926/FASTQ_2 \
                 --sample=SRR11537950 \
                 --localcores=56
