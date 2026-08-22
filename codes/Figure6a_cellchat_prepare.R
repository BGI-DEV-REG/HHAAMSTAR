# ============================================================
# Figure 6a - CellChat AGA Data Preparation
# Source: cellchat_group_submit-AGA.ipynb
# Output: data/Figure6a_cellchat.rds, Figure6a_cellchat_nep.csv, Figure6a_cellchat_net.csv
# ============================================================

options(stringsAsFactors = FALSE)
library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
library(viridis)
library(RColorBrewer)
library(ComplexHeatmap)

DATA_DIR <- "data"
OUT_DIR  <- "data"

obj <- readRDS(file.path(DATA_DIR, "HF_RNA_postQC_rpca_241011.rds"))
meta <- readRDS(file.path(DATA_DIR, "HF_RNA_rpca_meta_240612.rds"))
obj@meta.data <- meta[rownames(obj@meta.data), ]

obj <- subset(obj, sample %in% c("HB_Young", "HB_Mid_age", "AF", "AB"))

HB <- subset(obj, sample_type == "HB")
AB <- subset(obj, sample_type == "AB")
AF <- subset(obj, sample_type == "AF")

# 1. Create CellChat objects
cellchat_1 <- createCellChat(HB@assays$RNA@data, meta = HB@meta.data, group.by = "Abbreviation")
cellchat_2 <- createCellChat(AB@assays$RNA@data, meta = AB@meta.data, group.by = "Abbreviation")
cellchat_3 <- createCellChat(AF@assays$RNA@data, meta = AF@meta.data, group.by = "Abbreviation")

# 2. Set database
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat_1@DB <- CellChatDB.use
cellchat_2@DB <- CellChatDB.use
cellchat_3@DB <- CellChatDB.use

# 3. Preprocessing
cellchat_1 <- subsetData(cellchat_1)
cellchat_1 <- identifyOverExpressedGenes(cellchat_1)
cellchat_1 <- identifyOverExpressedInteractions(cellchat_1)
cellchat_1 <- projectData(cellchat_1, PPI.human)

cellchat_2 <- subsetData(cellchat_2)
cellchat_2 <- identifyOverExpressedGenes(cellchat_2)
cellchat_2 <- identifyOverExpressedInteractions(cellchat_2)
cellchat_2 <- projectData(cellchat_2, PPI.human)

cellchat_3 <- subsetData(cellchat_3)
cellchat_3 <- identifyOverExpressedGenes(cellchat_3)
cellchat_3 <- identifyOverExpressedInteractions(cellchat_3)
cellchat_3 <- projectData(cellchat_3, PPI.human)

# 4. Infer communication network
cellchat_1 <- computeCommunProb(cellchat_1, raw.use = FALSE, population.size = TRUE)
cellchat_1 <- filterCommunication(cellchat_1, min.cells = 10)
cellchat_1 <- computeCommunProbPathway(cellchat_1)
cellchat_1 <- aggregateNet(cellchat_1)
cellchat_1 <- netAnalysis_computeCentrality(cellchat_1, slot.name = "netP")

cellchat_2 <- computeCommunProb(cellchat_2, raw.use = FALSE, population.size = TRUE)
cellchat_2 <- filterCommunication(cellchat_2, min.cells = 10)
cellchat_2 <- computeCommunProbPathway(cellchat_2)
cellchat_2 <- aggregateNet(cellchat_2)
cellchat_2 <- netAnalysis_computeCentrality(cellchat_2, slot.name = "netP")

cellchat_3 <- computeCommunProb(cellchat_3, raw.use = FALSE, population.size = TRUE)
cellchat_3 <- filterCommunication(cellchat_3, min.cells = 10)
cellchat_3 <- computeCommunProbPathway(cellchat_3)
cellchat_3 <- aggregateNet(cellchat_3)
cellchat_3 <- netAnalysis_computeCentrality(cellchat_3, slot.name = "netP")

# 5. Merge and save
cco.list <- list(HB = cellchat_1, AB = cellchat_2, AF = cellchat_3)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)
saveRDS(cellchat, file.path(OUT_DIR, "Figure6a_cellchat.rds"))

# 6. Export communication tables
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
nep1 <- df.netp[["HB"]]; nep1$sample <- "HB"
nep2 <- df.netp[["AB"]]; nep2$sample <- "AB"
nep3 <- df.netp[["AF"]]; nep3$sample <- "AF"
nep <- rbind(nep1, nep2, nep3)
nep$p_log <- -log10(nep$pval)
nep$p_log <- ifelse(nep$p_log == Inf, 3, nep$p_log)
nep <- nep %>% mutate(source_target = paste(source, target, sep = "->"))
write.csv(nep, file.path(OUT_DIR, "Figure6a_cellchat_nep.csv"))

df.net <- subsetCommunication(cellchat, slot.name = "net")
net1 <- df.net[["HB"]]; net1$sample <- "HB"
net2 <- df.net[["AB"]]; net2$sample <- "AB"
net3 <- df.net[["AF"]]; net3$sample <- "AF"
net <- rbind(net1, net2, net3)
net$p_log <- -log10(net$pval)
net$p_log <- ifelse(net$p_log == Inf, 3, net$p_log)
net <- net %>% mutate(source_target = paste(source, target, sep = "->"))
write.csv(net, file.path(OUT_DIR, "Figure6a_cellchat_net.csv"))
