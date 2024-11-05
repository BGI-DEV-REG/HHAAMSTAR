rm(list=ls())
gc()

library(Seurat)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
setwd('../data/figure4/GO')

options(future.globals.maxSize = 300 * 1024^3)

ids = read.delim('/jdfssz1/ST_SUPERCELLS/P21Z10200N0171/Project/Hair_follicle/06.C4_scRNA/10.HFSC_AGAvsHB/HF_RNA_Spatial_map_hair.txt', header = F) %>%
    rename('id' = 'V1')
ids = ids$id

obj = lapply(ids, function(id){
    sample = gsub('_SCT_selection_\\d+','',id)
    obj = readRDS(paste0('../data/figure1/Spatial/04.RCTD/',sample,'/',id,'/',id,'.SCT.RCTD.rds'))
    obj$id_orig = id
    obj$id = gsub('selection_','',id)
    obj$sample = sample
    obj$cellid = colnames(obj)
    obj = RenameCells(obj, new.names = paste0(id, '__', Cells(obj)))
    obj@images = list()
    return(obj)
}) %>% Reduce(merge, .)
saveRDS(obj, file = 'fig4_GOmodule.rds')

hfdsc = read.csv('../data/figure4/GO/hHFDSC_GO.csv')
hfdsc = hfdsc[, c('Description','Hits')]
head(hfdsc,3)

Idents(obj) = 'RCTD'
obj1 = subset(obj,idents = 'hHFDSC') %>% SCTransform(., assay = 'Spatial')
obj1

for (n in 1:nrow(hfdsc)){
    gene = hfdsc[n, ]$Hits %>% strsplit(.,'\\|') %>% unlist
    gene = gene[gene %in% rownames(obj1)] %>% list
    Description = hfdsc[n, ]$Description %>% gsub(' ','_',.) %>% gsub(' - ','_',.) %>% gsub('-','_',.) %>% gsub('\n','_',.)
    if(length(gene[[1]]) < 1) next
    else{obj1 = AddModuleScore(obj1, features = gene, name = Description)}
}

qhfsc = read.csv('../data/figure4/GO/qHFSC_GO.csv')
qhfsc = qhfsc[, c('Description','Hits')]
head(qhfsc,3)

Idents(obj) = 'RCTD'
obj2 = subset(obj,idents = 'qHFSC') %>% SCTransform(., assay = 'Spatial')
obj2

for (n in 1:nrow(qhfsc)){
    gene = qhfsc[n, ]$Hits %>% strsplit(.,'\\|') %>% unlist
    gene = gene[gene %in% rownames(obj1)] %>% list
    Description = qhfsc[n, ]$Description %>% gsub(' ','_',.) %>% gsub(' - ','_',.) %>% gsub('-','_',.) %>% gsub('\n','_',.)
    if(length(gene[[1]]) < 1) next
    else{obj2 = AddModuleScore(obj2, features = gene, name = Description)}
}

for(go in colnames(obj1@meta.data)[16:length(colnames(obj1@meta.data))]){
    go1 = paste0('hfDSC_',go)
    obj@meta.data[, go1] = obj1@meta.data[Cells(obj), go]-min(obj1@meta.data[, go])
    obj@meta.data[, go1][is.na(obj@meta.data[, go1])] = 0
}

for(go in colnames(obj2@meta.data)[16:length(colnames(obj2@meta.data))]){
    go1 = paste0('qHFSC_',go)
    obj@meta.data[, go1] = obj2@meta.data[Cells(obj), go]-min(obj2@meta.data[, go])
    obj@meta.data[, go1][is.na(obj@meta.data[, go1])] = 0
}
saveRDS(obj@meta.data, file = 'fig4_module_meta.rds')

for (id in unique(obj$id_orig)){
    df = obj@meta.data[obj$id_orig == id, ]
    write.csv(df, file = paste0('../data/figure4/GO/Spatial_MarkerModule_meta/',id,'_fig4_module.csv'))
}

