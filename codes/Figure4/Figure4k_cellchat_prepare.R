# ============================================================
# Figure 4k - CellChat Age Data Preparation
# Source: cellchat_group_submit-Age.ipynb
# Output: data/Figure4k_cellchat.rds, Figure4k_cellchat_nep.csv, Figure4k_cellchat_net.csv
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

Idents(obj) <- "sample_type"
HB <- subset(obj, idents = "HB")

Idents(HB) <- "batch1"
batch_values <- c("5-HB-3y", "12-HB-55y", "16-HB-67y", "17-HB-66y")
HB <- subset(HB, idents = batch_values, invert = TRUE)

HB$Age <- as.numeric(as.character(HB$Age))
HB$group <- ifelse(HB$Age <= 30, "Young", ifelse(HB$Age >= 60, "Old", ""))

young <- subset(HB, subset = group == "Young")
old <- subset(HB, subset = group == "Old")

# 1. Create CellChat objects
cellchat_young <- createCellChat(young@assays$RNA@data, meta = young@meta.data, group.by = "Abbreviation")
cellchat_old <- createCellChat(old@assays$RNA@data, meta = old@meta.data, group.by = "Abbreviation")

# 2. Set database
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat_young@DB <- CellChatDB.use
cellchat_old@DB <- CellChatDB.use

# 3. Preprocessing
cellchat_young <- subsetData(cellchat_young)
cellchat_young <- identifyOverExpressedGenes(cellchat_young)
cellchat_young <- identifyOverExpressedInteractions(cellchat_young)
cellchat_young <- projectData(cellchat_young, PPI.human)

cellchat_old <- subsetData(cellchat_old)
cellchat_old <- identifyOverExpressedGenes(cellchat_old)
cellchat_old <- identifyOverExpressedInteractions(cellchat_old)
cellchat_old <- projectData(cellchat_old, PPI.human)

# 4. Infer communication network
cellchat_young <- computeCommunProb(cellchat_young, raw.use = FALSE, population.size = TRUE)
cellchat_young <- filterCommunication(cellchat_young, min.cells = 10)
cellchat_young <- computeCommunProbPathway(cellchat_young)
cellchat_young <- aggregateNet(cellchat_young)
cellchat_young <- netAnalysis_computeCentrality(cellchat_young, slot.name = "netP")

cellchat_old <- computeCommunProb(cellchat_old, raw.use = FALSE, population.size = TRUE)
cellchat_old <- filterCommunication(cellchat_old, min.cells = 10)
cellchat_old <- computeCommunProbPathway(cellchat_old)
cellchat_old <- aggregateNet(cellchat_old)
cellchat_old <- netAnalysis_computeCentrality(cellchat_old, slot.name = "netP")

# 5. Merge and save
cco.list <- list(Young = cellchat_young, Old = cellchat_old)
cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)
saveRDS(cellchat, file.path(OUT_DIR, "Figure4k_cellchat.rds"))

# 6. Export communication tables
df.net <- subsetCommunication(cellchat, slot.name = "net")
net1 <- df.net[["Old"]];   net1$age_group <- "Old"
net2 <- df.net[["Young"]]; net2$age_group <- "Young"
net <- rbind(net1, net2)
net$p_log <- -log10(net$pval)
net$p_log <- ifelse(net$p_log == Inf, 3, net$p_log)
net <- net %>% mutate(source_target = paste(source, target, sep = "->"))
write.csv(net, file.path(OUT_DIR, "Figure4k_cellchat_net.csv"))

df.netp <- subsetCommunication(cellchat, slot.name = "netP")
nep1 <- df.netp[["Old"]];   nep1$age_group <- "Old"
nep2 <- df.netp[["Young"]]; nep2$age_group <- "Young"
nep <- rbind(nep1, nep2)
nep$p_log <- -log10(nep$pval)
nep$p_log <- ifelse(nep$p_log == Inf, 3, nep$p_log)
nep <- nep %>% mutate(source_target = paste(source, target, sep = "->"))
write.csv(nep, file.path(OUT_DIR, "Figure4k_cellchat_nep.csv"))
