#load libraries

library(dplyr)
library(Seurat)
library(patchwork)

#load data
COVID_HUVEC.data <- Read10X(data.dir = '~//Documents/Data/COVID/Controls_HUVECs/primary_alignment_all_wCtrl/outs/filtered_feature_bc_matrix')

# create Seurat object
COVID_HUVEC <- CreateSeuratObject(counts = COVID_HUVEC.data, project = 'COVID_HUVEC', min.cells = 3, min.features = 200)
COVID_HUVEC

# inspect metadata
head(COVID_HUVEC@meta.data)
tail(COVID_HUVEC@meta.data)


# add barcode metadata - extract patient ID from barcode by strpsplit on '-' and extracting second element of resulting list (here the indentifier is appended after '-' to the barcode)
Sample_ID <- sapply(strsplit(rownames(COVID_HUVEC@meta.data), split = '-'), "[[",2)

# add barcode metadata - supply pateint ID as additional metadata
COVID_HUVEC <- AddMetaData(object = COVID_HUVEC, metadata = data.frame(Sample_ID = Sample_ID, row.names = rownames(COVID_HUVEC@meta.data)))

#inspect
head(COVID_HUVEC@meta.data)
tail(COVID_HUVEC@meta.data)

COVID_HUVEC@meta.data$Sample_ID = dplyr::recode(COVID_HUVEC@meta.data$Sample_ID,
																								"1"="OC43_12H",
																								"2"="OC43_24H",
																								"3" ="CTRL")
#QC 
COVID_HUVEC[['percent.mt']] <- PercentageFeatureSet(COVID_HUVEC, pattern = '^MT-')

Idents(COVID_HUVEC) <- 'Sample_ID'
VlnPlot(COVID_HUVEC, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol=3)

COVID_HUVEC_QC <- subset(COVID_HUVEC, subset =  nCount_RNA < 50000 & percent.mt < 15)
VlnPlot(COVID_HUVEC_QC, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol=3)

#normalization etc
COVID_HUVEC <- NormalizeData(COVID_HUVEC)
COVID_HUVEC <- FindVariableFeatures(COVID_HUVEC)
top10 <- head(VariableFeatures(COVID_HUVEC), 10)
plot1 <- VariableFeaturePlot(COVID_HUVEC)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

all.genes <- rownames(COVID_HUVEC)
COVID_HUVEC <- ScaleData(COVID_HUVEC, features = all.genes)
COVID_HUVEC <- RunPCA(COVID_HUVEC)
ElbowPlot(COVID_HUVEC)
DimHeatmap(COVID_HUVEC, dims= 1:20, cells = 500, balanced = TRUE)

COVID_HUVEC <- RunUMAP(COVID_HUVEC, dims = 1:20, reduction = 'pca')
COVID_HUVEC <- FindNeighbors(COVID_HUVEC, dims = 1:20)

DimPlot(COVID_HUVEC, group.by = 'Sample_ID')

# hCoV-OC43 (NC-006213-1) expression
FeaturePlot(COVID_HUVEC, features = 'NC-006213-1')

VlnPlot(COVID_HUVEC, features = 'NC-006213-1', group.by = 'Sample_ID', sort = 'decreasing')


