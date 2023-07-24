library(CellChat)
library(patchwork)
library(Seurat)
library(NMF)
library(ggalluvial)
library(future)


options(stringsAsFactors = FALSE)
options(future.globals.maxSize= 6291456000)

args <- commandArgs(trailingOnly = TRUE)

obj <- readRDS(args[1])
num_clusters <- length(unique(sort(obj@meta.data$region)))
data.input <- as.matrix(obj@assays$Spatial@data)
meta_data <- cbind(rownames(obj@meta.data), obj@meta.data[,'region', drop=F])
meta_data <- as.matrix(meta_data)

cellchat <- createCellChat(object = data.input, meta = meta_data, group.by = "region")

if ( args[2] == 'mouse' ) {
    CellChatDB <- CellChatDB.mouse
} else if ( args[2] == 'human' ) {
    CellChatDB <- CellChatDB.human
}
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat)
## subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#cellchat <- projectData(cellchat, PPI.mouse)   ##optional

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net,paste(args[3],"net_legand_receptor.csv",sep=""))

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name ="netP")
write.csv(df.netp,paste(args[3],"net_pathway.csv",sep=""))


cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

pdf('cellchat_network.pdf', width = 9, height = 9)
#par(mfrow = c(1,2))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
#netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

png('cellchat_network.png', width = 600, height = 600)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

mat <- cellchat@net$weight
pdf('cellchat_seperate.pdf', width = 20, height = ceiling(num_clusters/4)*5 )
par(mfrow = c(ceiling(num_clusters/4),4))
for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()
png('cellchat_seperate.png', width = 800, height = ceiling(num_clusters/4)*200 )
par(mfrow = c(ceiling(num_clusters/4),4))
for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

saveRDS(cellchat,paste(args[3],".cellchat.rds",sep=""))

