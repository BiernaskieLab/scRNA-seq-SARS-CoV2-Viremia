
#! /usr/bin/env bash

#BSUB -J CR3.0.0
#BSUB -n 56
#BSUB -R "span[hosts=1]"
#BSUB -W 599:59
#BSUB -o COVID_72_Hour_Cov3%J.out
#BSUB -e COVID_72_Hour_Cov3%J.err
#BSUB -N

/home/ajaffer/cellranger_3.0.1/cellranger-3.0.1/cellranger aggr --id=Agg_COVID_72_Hour_Cov3 \
                --csv=Agg.csv \
                --normalize=mapped
