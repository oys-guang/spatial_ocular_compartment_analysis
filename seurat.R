library(Seurat)
library(ggplot2)
library(patchwork)
library(plyr)
library(dplyr)
library(tidyr)
library(data.table)
library(Matrix.utils)
library(future)
library(data.table)
library(reticulate)
args <- commandArgs(trailingOnly = TRUE)
##########
#read in data
t1<-Sys.time()

infile=args[1]

data<-readRDS(infile)
obj<-subset(data, subset = nCount_Spatial >= 20 & nFeature_Spatial >= 20)
#plot1 <- VlnPlot(data, features =c("nCount_Spatial","nFeature_Spatial"), pt.size = 0.1) + NoLegend()
#plot2 <- SpatialFeaturePlot(obj, stroke = 0,pt.size.factor = 10,features = c("nCount_Spatial","nFeature_Spatial")) + theme(legend.position = "right")
#pdf('vln.pdf', width = 8, height = 9)
#plot1
#dev.off()
#pdf('expr.pdf', width = 8, height = 9)
#plot2
#dev.off()
obj <- SCTransform(obj, assay = "Spatial", verbose = FALSE)

#path_to_python = "/dellfsqd2/ST_OCEAN/USER/maoxuebin/5.software/anaconda3/envs/conda/bin/python"
#use_python(path_to_python)
#genes <- as.data.frame(rownames(obj), row.names = rownames(obj))
#names(genes) <- "Gene"

#cells <- as.data.frame(colnames(obj), row.names = colnames(obj))
#names(cells) <- "CellID"

#row <- obj@images[[1]]@coordinates$row
#col <- obj@images[[1]]@coordinates$col
#coordinates <- list(matrix(c(row, col), ncol = 2))
#names(coordinates) <- "spatial"

#ann <- import("anndata")
#ad <- ann$AnnData(X = obj@assays$SCT@data, obs = genes, var = cells, varm = coordinates,layers = list(counts = obj@assays$SCT@counts))
#ad <- ad$T

#ad$write_h5ad("12L_SCT_raw.h5ad")


obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)
#pdf("PCA.pdf",width=12,height=9)
#ElbowPlot(obj)
#dev.off
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:20)
obj <- FindClusters(obj, verbose = FALSE,resolution =0.5)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:20)
p1 <- DimPlot(obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(obj,stroke=0,pt.size.factor=20)
pdf('umap.pdf', width = 12, height = 9)
p1
dev.off()
pdf('spatial.pdf', width = 30, height = 30)
p2
dev.off()
saveRDS(obj,file="seurat.rds")

