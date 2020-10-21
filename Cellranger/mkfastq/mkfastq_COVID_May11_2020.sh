
#! /usr/bin/env bash

#BSUB -J CR3.0.0
#BSUB -n 56
#BSUB -R "span[hosts=1]"
#BSUB -W 999:59
#BSUB -o COVID_May_11%J.out
#BSUB -e COVID_May_11%J.err
#BSUB -N




export PATH=/export/common/programs/bcl2fastq2-v2.20/bin:$PATH
/home/ajaffer/cellranger_3.0.1/cellranger-3.0.1/cellranger mkfastq --id=1_COVID_B3 \
                     --run=/home/ajaffer/NovaSeq_May2020 \
                     --localcores=56 \
                     --jobmode=local \
                     --simple-csv=/home/ajaffer/COVID_May11_2020/cellranger-mkfastq-COVID_1_Blinded_B3.csv



/home/ajaffer/cellranger_3.0.1/cellranger-3.0.1/cellranger mkfastq --id=2_COVID_B4 \
                     --run=/home/ajaffer/NovaSeq_May2020 \
                     --localcores=56 \
                     --jobmode=local \
                     --simple-csv=/home/ajaffer/COVID_May11_2020/cellranger-mkfastq-COVID_2_Blinded_B4.csv



/home/ajaffer/cellranger_3.0.1/cellranger-3.0.1/cellranger mkfastq --id=3_COVID_B5 \
                     --run=/home/ajaffer/NovaSeq_May2020 \
                     --localcores=56 \
                     --jobmode=local \
                     --simple-csv=/home/ajaffer/COVID_May11_2020/cellranger-mkfastq-COVID_3_Blinded_B5.csv



/home/ajaffer/cellranger_3.0.1/cellranger-3.0.1/cellranger mkfastq --id=4_COVID_B6 \
                     --run=/home/ajaffer/NovaSeq_May2020 \
                     --localcores=56 \
                     --jobmode=local \
                     --simple-csv=/home/ajaffer/COVID_May11_2020/cellranger-mkfastq-COVID_4_Blinded_B6.csv



/home/ajaffer/cellranger_3.0.1/cellranger-3.0.1/cellranger mkfastq --id=5_COVID_B7 \
                     --run=/home/ajaffer/NovaSeq_May2020 \
                     --localcores=56 \
                     --jobmode=local \
                     --simple-csv=/home/ajaffer/COVID_May11_2020/cellranger-mkfastq-COVID_5_Blinded_B7.csv



/home/ajaffer/cellranger_3.0.1/cellranger-3.0.1/cellranger mkfastq --id=6_COVID_B8 \
                     --run=/home/ajaffer/NovaSeq_May2020 \
                     --localcores=56 \
                     --jobmode=local \
                     --simple-csv=/home/ajaffer/COVID_May11_2020/cellranger-mkfastq-COVID_6_Blinded_B8.csv




/home/ajaffer/cellranger_3.0.1/cellranger-3.0.1/cellranger mkfastq --id=7_COVID_B9 \
                     --run=/home/ajaffer/NovaSeq_May2020 \
                     --localcores=56 \
                     --jobmode=local \
                     --simple-csv=/home/ajaffer/COVID_May11_2020/cellranger-mkfastq-COVID_7_Blinded_B9.csv



/home/ajaffer/cellranger_3.0.1/cellranger-3.0.1/cellranger mkfastq --id=8_COVID_B10 \
                     --run=/home/ajaffer/NovaSeq_May2020 \
                     --localcores=56 \
                     --jobmode=local \
                     --simple-csv=/home/ajaffer/COVID_May11_2020/cellranger-mkfastq-COVID_8_Blinded_B10.csv






