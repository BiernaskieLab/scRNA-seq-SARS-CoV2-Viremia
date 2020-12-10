#load libraries

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)

#load data
COVID1.data <- Read10X(data.dir = '~/desktop/COVID_Aggregated_Files_2/Aggr_72_Hour_Cov2_BALF/filtered_feature_bc_matrix')

# create Seurat object
COVID1 <- CreateSeuratObject(counts = COVID1.data, project = 'COVID1', min.cells = 3, min.features = 200)
COVID1

# inspect metadata
head(COVID1@meta.data)
tail(COVID1@meta.data)


# add barcode metadata - extract patient ID from barcode by strpsplit on '-' and extracting second element of resulting list (here the indentifier is appended after '-' to the barcode)
Sample_ID <- sapply(strsplit(rownames(COVID1@meta.data), split = '-'), "[[",2)

# add barcode metadata - supply pateint ID as additional metadata
COVID1 <- AddMetaData(object = COVID1, metadata = data.frame(Sample_ID = Sample_ID, row.names = rownames(COVID1@meta.data)))

#inspect
head(COVID1@meta.data)
tail(COVID1@meta.data)

#change names of Sample_ID

COVID1@meta.data$Sample_ID = dplyr::recode(COVID1@meta.data$Sample_ID,
																										"1"="UC_1",
																										"2"="UC_4",
																										"3"="UC_5",
																										"4"="UC_6",
																										"5"="UC_7",
																										"6"="UC_8",
																										"7"="BALF_54_mild",
																										"8"="BALF_55_mild",
																										"9"="BALF_56_severe",
																										"10" ="BALF_57_mild",
																										"11" ="BALF_58_severe",
																										"12" ="BALF_59_severe",
																										"13" ="BALF_46_CD45_HC",
																										"14" ="BALF_47_CD45_HC",
																										"15" ="BALF_48_CD45_HC",
																										"16" ="BALF_49_severe",
																										"17" ="BALF_50_severe",
																										"18" ="BALF_51_severe")


#subset to remove UC_1 (low cell number and sub-par collection)

Idents(COVID1) <- 'Sample_ID'

COVID1_subset <- SubsetData(COVID1, ident.remove = 'UC_1' )
head(COVID1_subset@meta.data)

#QC 
COVID1_subset[['percent.mt']] <- PercentageFeatureSet(COVID1_subset, pattern = '^MT-')

VlnPlot(COVID1_subset, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol=3)

COVID1_QC <- subset(COVID1_subset, subset = nFeature_RNA < 5000 & nCount_RNA < 30000 & percent.mt < 10)
VlnPlot(COVID1_QC, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol=3)


#normalization etc
COVID1 <- NormalizeData(COVID1_QC)
COVID1 <- FindVariableFeatures(COVID1)
top10 <- head(VariableFeatures(COVID1), 10)
plot1 <- VariableFeaturePlot(COVID1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# scale data - run with features = all.genes if sufficient compute power 
# all.genes <- rownames(COVID1)
# COVID1 <- ScaleData(COVID1, features = all.genes)
COVID1 <- ScaleData(COVID1)
COVID1 <- RunPCA(COVID1)
ElbowPlot(COVID1)
COVID1 <- RunUMAP(COVID1, dims = 1:20, reduction = 'pca')
COVID1 <- FindNeighbors(COVID1, dims = 1:20)
COVID1 <- FindClusters(COVID1, resolution = 1.0)

# analysis
DimPlot(COVID1, group.by = 'Sample_ID', label = FALSE, order = 'decreasing')
DimPlot(COVID1, group.by = 'seurat_clusters', label = TRUE)

head(COVID1@meta.data)

# SARS-CoV-2 detection, receptor expression by Sample_ID
VlnPlot(COVID1, features = 'SARS-CoV-2', group.by = 'Sample_ID', sort = 'decreasing')
VlnPlot(COVID1, features = 'ACE2', group.by = 'Sample_ID', sort = 'decreasing')
VlnPlot(COVID1, features = 'TMPRSS2', group.by = 'Sample_ID', sort = 'decreasing')
VlnPlot(COVID1, features = 'NRP1', group.by = 'Sample_ID', sort = 'decreasing')

FeaturePlot(COVID1, features = 'SARS-CoV-2', max.cutoff = 4)
FeaturePlot(COVID1, features = c('ACE2','TMPRSS2','NRP1','MRC1'))

# cluster annotation
COVID1.markers <- FindAllMarkers(COVID1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(COVID1.markers,'COVID128cluster_findallmarkers.csv')

top10clustermarkers <- COVID1.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
write.csv(top10clustermarkers, 'COVID1_top10markers.csv')

FeaturePlot(COVID1, features = c('CD14', 'LYZ', 'CD4', 'CD8A', 'CD68', 'IL7R', 'CCR7', 'S100A4','NKG7','PPBP','CST3','FCER1A','MS4A1', 'KRT14','MPO','FCGR3A'))
FeaturePlot(COVID1, features = c('CD19', 'CD3G', 'TRGC2', 'TRAC',  'TRBC2',  'CD3D',  'CD27', 'FOXP3','IL2RA','PTPRC', 'ITGAM','ITGAX','IGHD','CD38', 'CD74', 'IGHD', 'NCAM1'))
# broad
FeaturePlot(COVID1, features = c('CD19', 'CD3G','IL2RA','PTPRC', 'ITGAM','ITGAX','IGHD','CD38', 'CD74', 'IGHD', 'NCAM1'))

# neutrophils - from Elodie
FeaturePlot(COVID1, features = c('CXCR2', 'TLR4', 'MMP9', 'CD177', 'FCGR3B', 'CSF3R', 'S100A8', 'IL1R2', 'ALPL', 'ANXA3', 'PROK2'))
# B cells
FeaturePlot(COVID1, features = c('CD19','IGHD', 'CD34','PTPRC','CD27', 'CD38'))
# NK cells
FeaturePlot(COVID1, features = c('NCAM1', 'NKG7', 'GNLY', 'FCGR3A', 'FCGR3B'))
# T cell subtypes - CD4 CD8 TReg Th17
FeaturePlot(COVID1, features = c('CD3G', 'TRGC2', 'TRAC',  'TRBC2',  'CD3D', 'FOXP3','IL2RA','PTPRC', 'CD74', 'CD4', 'CD8A'))
# CD4 Tcells - from Szabo_2019 nature communications paper
FeaturePlot(COVID1, features = c('CCR7','SELL', 'TCF7', 'IL2', 'IFNG', 'IL4', 'IL13', 'IL17A', 'TNF', 'PRF1'))
# monocytes
FeaturePlot(COVID1, features = c('ITGAM','ITGAX','CD14', 'FCGR3A', 'FCGR3B'))
# macrophage subsets
FeaturePlot(COVID1, features = c('CD14','FCGR3A', 'FCGR3B', 'PTPRC', 'ITGAM', 'CD68', 'MRC1', 'ITGAX', 'SIGLEC1'))
# DCs
FeaturePlot(COVID1, features = c('PTPRC', 'CD14','IL3RA','ITGAX', 'CD74'))
# airway epithelial cells
FeaturePlot(COVID1, features = c('KRT4','EPCAM', 'CDH1','AGER','HOPX','SFTPC', 'SLC34A2', 'ABCA3', 'SOX2', 'PAX9','TP63','KRT5', 'MUC5B','SCGB1A1'))
# dividing 
FeaturePlot(COVID1, features = c('MKI67'))



# VlnPlots 
VlnPlot(COVID1, features = 'FOXP3', group.by = 'seurat_clusters')
VlnPlot(COVID1, features = 'CD19', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'IGHD', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'CD4', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'NCAM1', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'CD74', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'ITGAM', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'ITGAX', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'FCGR3B', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'FCGR3A', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'CCR3', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'MBP', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'IL2RA', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'CD19', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'IGLV3-19', group.by = 'seurat_clusters', sort = 'decreasing')
VlnPlot(COVID1, features = 'PDGFRA', group.by = 'seurat_clusters', sort = 'decreasing')

# add barcode metadata - cell type annotation - COV2
new.cluster.ids <- c("Mac_1","Mac_1", "Mac_1", "Mac_2","UC_TCells",
										 "Mac_2", "Neutrophils", "TCells", "Mac_2", "Neutrophils", 
										 "Mac_1", "UC_TCells", "Mac_1", "UC_TCells", "TCells",
										 "Neutrophils", "UC_BCells", "Mac_2", "Mac_1", "UC_TCells", 
										 "TCells", "BCells", "DC", "TCells", "Epithelial", 
										 "NK", "Epithelial", "NK", "UC_Monocytes", "UC_TCells",
										 "Mac_1", "UC_Platelets","BCells", "Mac_2", "Mac_2",
										 "Epithelial", "Myeloid", "BCells")

Idents(COVID1) <- 'seurat_clusters'
names(new.cluster.ids) <- levels(COVID1)
COVID1 <- RenameIdents(COVID1, new.cluster.ids)
head(COVID1@meta.data)
head(Idents(COVID1))

# save new idents in metadata
COVID1[['short_clusterID']] <- Idents(COVID1)

Idents(COVID1) <- 'short_clusterID'

DimPlot(COVID1, group.by = 'short_clusterID', order = c('UC_BCells', 'UC_TCells', 'UC_Monocytes', 'UC_Platelets'), label = FALSE)

VlnPlot(COVID1, features = 'SARS-CoV-2', group.by = 'short_clusterID', sort = 'decreasing')
VlnPlot(COVID1, features = 'ACE2', group.by = 'short_clusterID', sort = 'decreasing')
VlnPlot(COVID1, features = 'TMPRSS2', group.by = 'short_clusterID', sort = 'decreasing')
VlnPlot(COVID1, features = 'NRP1', group.by = 'short_clusterID', sort = 'decreasing')

Idents(COVID1) <- 'Sample_ID'
table(Idents(COVID1))

table(COVID1@meta.data$short_clusterID, COVID1@meta.data$CoV2_status)

Expression <- AverageExpression(COVID1, features = c('SARS-CoV-2', 'ACE2', 'TMPRSS2', 'NRP1'), verbose = TRUE)
head(Expression)
write.csv(Expression, 'COVID1_4geneExpression.csv')


# subset SARS-CoV-2 +ve
COVID1_CoV2 <- subset(COVID1, subset = `SARS-CoV-2` > 0.1)
table(Idents(COVID1_CoV2))
saveRDS(COVID1_CoV2, 'V2Cov2_COVID1_CoV2.rds')

# new ident for CoV2 +ve

COVID1 <- SetIdent(object = COVID1,  value = 'CoVNeg')
COVID1 <- SetIdent(object = COVID1, cells = Cells(COVID1_CoV2), value = 'CoVPos')

table(Idents(COVID1))
CoV2_status <- c('CoVPos','CoVNeg')
names(CoV2_status) <- levels(COVID1)
COVID1 <- RenameIdents(COVID1, CoV2_status)
COVID1[['CoV2_status']] <- Idents(COVID1)
head(COVID1@meta.data)
tail(COVID1@meta.data)

saveRDS(COVID1, 'V2Cov2_COVID1_CoVstatus.rds')

vlnplot(COVID1, features = 'SARS-CoV-2', group.by = 'short_clusterID', split.by = 'CoV2_status')
