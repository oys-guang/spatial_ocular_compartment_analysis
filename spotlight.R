library(Matrix)
library(data.table)
library(Seurat)
library(SeuratData)
library(dplyr)
library(gt)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
library(png)
library(ggsci)
args <- commandArgs(trailingOnly = TRUE)

st <- readRDS(args[1])
img <- GetImage(st, mode='raw')
writePNG(img, target='log.png')

sc <- readRDS(args[2])

set.seed(123)
sc <- Seurat::SCTransform(sc, verbose = FALSE) %>% Seurat::RunPCA(., verbose = FALSE) %>% Seurat::RunUMAP(., dims = 1:30, verbose = FALSE)
Seurat::Idents(object = sc) <- sc@meta.data$celltype
cluster_markers_all <- Seurat::FindAllMarkers(object = sc, 
                                              assay = "SCT",
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE,
					      min.pct = 0.25,
                                              logfc.threshold = 0.25
                                              )

topn <- cluster_markers_all %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

spotlight_ls <- spotlight_deconvolution(
                    se_sc = sc,
                    counts_spatial = st@assays$Spatial@counts,
                    clust_vr = "celltype", # Variable in sc containing the cell-type annotation
#                    cluster_markers = cluster_markers_all, # Dataframe with the marker genes
                    cluster_markers = topn,
                    cl_n = 5000, # number of cells per cell type to use
                    hvg = 3000, # Number of HVG to use
                    ntop = NULL, # How many of the marker genes to use (by default all)
                    transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
                    method = "nsNMF", # Factorization method
                    min_cont = 0.05 # Remove those cells contributing to a spot below a certain threshold 
                    )

nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

h <- NMF::coef(nmf_mod[[1]])
rownames(h) <- paste("Topic", 1:nrow(h), sep = "_")
topic_profile_plts <- SPOTlight::dot_plot_profiles_fun(
                                                       h = h,
                                                       train_cell_clust = nmf_mod[[2]]
                                                       )

topic_profile_plts[[2]] + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45), 
                                         axis.text = ggplot2::element_text(size = 16))
ggsave('./topic_profile.pdf', width = 240, height = 240, units = 'mm')

topic_profile_plts[[1]] + theme(axis.text.x = element_text(angle = 45), 
                                axis.text = element_text(size = 16))
ggsave('./topic_proportion.pdf', width = 720, height = 600, units = 'mm')

basis_spotlight <- data.frame(NMF::basis(nmf_mod[[1]]))

colnames(basis_spotlight) <- unique(stringr::str_wrap(nmf_mod[[2]], width = 30))

#basis_spotlight %>% dplyr::arrange(desc(Astro)) %>% round(., 5) %>% DT::datatable(., filter = "top")

# This is the equivalent to setting min_cont to 0.04
decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
rownames(decon_mtrx) <- colnames(st)

decon_df <- decon_mtrx %>% data.frame() %>% 
            tibble::rownames_to_column("barcodes")

st@meta.data <- st@meta.data %>% 
                tibble::rownames_to_column("barcodes") %>% 
                dplyr::left_join(decon_df, by = "barcodes") %>%
                tibble::column_to_rownames("barcodes")

#Seurat::SpatialFeaturePlot(object = st,
#                           features = c("EVM%20prediction%20000000F_arrow_pilon.106"),
#                           alpha = c(0.1, 1))

#saveRDS(decon_mtrx, 'spotlight.rds')

cell_types_all <- colnames(decon_mtrx)[which(colnames(decon_mtrx) != "res_ss")]
#cell_types_all <- paste('X', cell_types_all, sep='')

predict_CellType = apply(decon_mtrx[,which(colnames(decon_mtrx) != "res_ss")], 1, function(x){
                         index = which.max(x)
                         a = colnames(decon_mtrx)[index]
                         return(a)})

st@meta.data$predicted_CellType <- predict_CellType

saveRDS(st, 'spotlight.rds')

cluster_Palette <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')

meta.data <- st@meta.data
for (sample in unique(meta.data$orig.ident)){
    meta.data$x[meta.data$orig.ident == sample] <- meta.data$x[meta.data$orig.ident == sample] - min(meta.data$x[meta.data$orig.ident == sample])
    meta.data$y[meta.data$orig.ident == sample] <- meta.data$y[meta.data$orig.ident == sample] - min(meta.data$y[meta.data$orig.ident == sample])
}

meta.data %>% ggplot(aes_string(x = 'x', y = 'y', color = 'predicted_CellType')) +
        geom_point(shape = 19, size = 0.5) +
        facet_wrap(~ orig.ident, nrow = 2) +
        theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(),
                axis.title = element_blank(), axis.line = element_blank(), legend.position = 'right',
                panel.background = element_blank()) +
        scale_color_manual(values = cluster_Palette) +
        guides(colour = guide_legend(override.aes = list(size=8), nrow = 1)) +
        coord_fixed(ratio=1)

ggsave('./projection_feature.pdf', width = 100, height = 100, units = 'cm')

decon_mtrx_sub <- decon_mtrx[, cell_types_all]
decon_mtrx_sub <- decon_mtrx_sub[, colSums(decon_mtrx_sub) > 0]
decon_mtrx_cor <- cor(decon_mtrx_sub)
p.mat <- corrplot::cor.mtest(mat = decon_mtrx_sub, conf.level = 0.95)
ggcorrplot::ggcorrplot(
  corr = decon_mtrx_cor,
  p.mat = p.mat[[1]],
  hc.order = TRUE,
  type = "full",
  insig = "blank",
  lab = TRUE,
  outline.col = "lightgrey",
  method = "square",
  # colors = c("#4477AA", "white", "#BB4444"))
  colors = c("#6D9EC1", "white", "#E46726"),
  title = "Predicted cell-cell proportion correlation",
  legend.title = "Correlation\n(Pearson)") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(size = 22, hjust = 0.5, face = "bold"),
    legend.text = ggplot2::element_text(size = 12),
    legend.title = ggplot2::element_text(size = 15),
    axis.text.x = ggplot2::element_text(angle = 90),
    axis.text = ggplot2::element_text(size = 18, vjust = 0.5))
ggsave('./correlation.pdf', width = 400, height = 400, units = 'mm')
#SPOTlight::spatial_scatterpie(se_obj = st, 
#                              cell_types_all = cell_types_all, 
#                              img_path = "log.png", 
#                              pie_scale = 0.1)


