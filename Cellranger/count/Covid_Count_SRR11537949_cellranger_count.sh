#! /usr/bin/env bash

#BSUB -J CR3.0.0
#BSUB -n 56
#BSUB -R "span[hosts=1]"
#BSUB -W 999:59
#BSUB -o Covid_Count_SRR11537949%J.out
#BSUB -e Covid_Count_SRR11537949%J.err
#BSUB -N


/home/ajaffer/cellranger_3.0.1/cellranger-3.0.1/cellranger count --id=SRR11537949_May2020 \
                 --transcriptome=/home/ajaffer/refdata-cellranger-GRCh38-3.0.0_nCov_19/GRCh38_and_SARS-Cov2 \
                 --fastqs=/home/ajaffer/COVID_GEO_GSE145926/FASTQ_2 \
                 --sample=SRR11537949 \
                 --localcores=56

/home/ajaffer/cellranger_3.0.1/cellranger-3.0.1/cellranger count --id=SRR11537949_GRCh38only \
                 --transcriptome=/home/ajaffer/refdata-cellranger-GRCh38-3.0.0_nCov_19/refdata-cellranger-GRCh38-3.0.0 \
                 --fastqs=/home/ajaffer/COVID_GEO_GSE145926/FASTQ_2 \
                 --sample=SRR11537949 \
                 --localcores=56

