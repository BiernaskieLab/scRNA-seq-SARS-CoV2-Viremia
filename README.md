# scRNA-seq-SARS-CoV2-Viremia
This repository contains analysis scripts for Rosin et al. 2020 (Under Review).

# Summary

## Abstract
In late 2019 a novel coronavirus (SARS-CoV-2) emerged, and has since caused a global pandemic. Understanding the pathogenesis of COVID-19 disease is necessary to inform development of therapeutics, and management of infected patients. Using scRNAseq of blood drawn from SARS-CoV-2 patients, we asked whether SARS-CoV-2 may exploit immune cells as a ‘Trojan Horse’ to disseminate and access multiple organ systems. Our data suggests that circulating cells are not actively infected with SARS-CoV-2, and do not appear to be a source of viral dissemination.

# Data

## Single-cell RNA-Seq
NCBI GEO: [GSE151969](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151969) [GSE156639](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156639)<br/>
```
GSE151969:
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151969/suppl/GSE151969_RAW.tar
tar -xvf GSE151969_RAW.tar
GSE156639:
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156639/suppl/GSE156639_RAW.tar
tar -xvf GSE156639_RAW.tar
```
NCBI SRA: To be released <br/>
```
source activate sratoolkit
fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR_ID
## Specify SRR_ID - obtained using SRA Run selector.
```

# Toolkits used
`Cellranger v.3.1.0` - Alignment and aggregation of 10x-generated scRNA-Seq data.  <br/>
`Seurat v.3.1.5` - scRNA-Seq Analysis. <br/>
`R v.3.6.1 ` - scRNA-Seq Analysis. <br/>

# Contact
Dr. Jeff Biernaskie (jabierna@ucalgary.ca)<br/>
Dr. Nicole L. Rosin (Nicole.Rosin@ucalgary.ca)<br/>
Arzina Jaffer (arzina.jaffer1@ucalgary.ca)<br/>
Sarthak Sinha (sarthak.chinoo@gmail.com)
