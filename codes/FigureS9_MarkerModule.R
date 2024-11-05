rm(list=ls())
gc()

library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
setwd('../data/figure1/Spatial/05.MarkerModule')

options(future.globals.maxSize = 200 * 1024^3)

sc = readRDS('../data/figure1/scRNA-seq/HF_RNA_postQC_rpca.rds')
Idents(sc) = 'celltype'
marker = FindAllMarkers(sc, only.pos = T, logfc.threshold = 0.5)
write.csv(marker, file = '../data/figure1/scRNA-seq/HF_RNA_CellMarker.csv')

marker = read.csv('../data/figure1/scRNA-seq/HF_RNA_CellMarker.csv',row.names=1) %>% 
    filter(p_val < 0.05) %>% #filter(avg_log2FC > 1) %>%
    group_by(cluster2) %>%
    arrange(desc(avg_log2FC),p_val) %>%
    slice(1:10)

# before cell markers' module score calculate
# we use cirro to lasso each hair follicle of sclap tissue of Stereo-seq
# the object was subset from the whole tissue by cellid
    
sample = gsub('_SCT_selection_\\d+','',id)
obj = readRDS(paste0('../data/figure1/Spatial/04.RCTD/',sample,'/',id,'/',id,'.SCT.RCTD.rds'))
obj$id = id
obj$sample = sample
# obj = RenameCells(obj, new.names = paste0(id, '__', Cells(obj)))
DefaultAssay(obj) = 'SCT'

gene = marker[marker$cluster == celltype, ] %>% pull(gene)
gene = gene[gene %in% rownames(obj)] %>% list
obj = AddModuleScore(obj, features = gene, name = paste0(cluster,'_Features'),seed = 2024)
write.csv(obj@meta.data, file = paste0('../data/figure1/Spatial/05.MarkerModule/',sample,'/',id,'_',celltype,'_Spatial_MarkerModule.csv'))
