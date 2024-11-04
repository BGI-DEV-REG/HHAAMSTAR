library(Seurat)
library(ggplot2)
library(dplyr)
library(stringr)
library(grid)
library(tibble)
library(stringr)
library(parallel)

source('../data/figure1/code/seurat_helper.r')

# removing contamination RNAs

data_path='../data/figure1/scRNA-seq/raw_data/'
output_path='../data/figure1/scRNA-seq/SoupX/'
data_sample = read.csv('../data/figure1/scRNA-seq/Figure1B_scRNA-raw_data_info.csv',row.names=1)
data_names = data_sample$data_names
sample_names = data_sample$sample_names

# add matrix path
matrix_paths <- c()
matrix_path1 <- c()
matrix_path2 <- c()
output_paths <- c()
for (i in 1:length(data_names)){
    #
    matrix_paths=c(matrix_paths, paste0(data_path, data_names[i]))
    matrix_path1=c(matrix_path1, paste0(data_path, data_names[i], '/04.Matrix/FilterMatrix/'))
    matrix_path2=c(matrix_path2, paste0(data_path, data_names[i], '/02.cDNAAnno/RawMatrix/'))
    #
    output_paths=c(output_paths, paste0(output_path, data_names[i]))
}
# output data
for (i in 1:length(matrix_paths)){
    out <- runSoupX_Droplets(matrix_path=matrix_paths[i])
    #if(!dir.exists(output_paths[i])){dir.create(output_paths[i])}
    ExportData_10X_format(out,out_dir=output_paths[i], assay='RNA')
}

# QC

## add matrix path
matrix_paths <- c()
for (i in 1:length(data_names)){
    matrix_paths=c(matrix_paths, paste0(data_path, data_names[i]))
}

## merge samples 
### create seurat obj
seurat.list=c()
for (i in 1:length(matrix_paths)){
    seurat.obj=Create_scRNA_object(data.dir = matrix_paths[i])
    #
    seurat.obj$batch=sample_names[i]
    seurat.obj$tissue='scalp'
    seurat.obj$orig.ident=sample_names[i]
    #
    seurat.list=c(seurat.list, seurat.obj)
}

### get overlap genes
genes=Overlap_Seurat_Genes(seurat.list)
### extract genes of each seurat object
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
        x <- x[genes, ]
})

## QC
seurat.list <- lapply(X = seurat.list, FUN = function(x) {
    x <- QC_RNA(proj = x, genome = "hg38")
    # These defaults can be run just by providing accepted species name
    x <- Add_Mito_Ribo_Seurat(seurat_object = x, species = "Human")
    # These defaults can be run just by providing accepted species name
})

## outliers and doublets
for (i in 1:length(sample_names)){
    # outliers
    seurat.list[[i]] <- subset(seurat.list[[i]], subset = nCount_RNA > 100)
    outlier <- is_outlier(seurat.list[[i]], "nCount_RNA", 4) | 
                is_outlier(seurat.list[[i]], "nFeature_RNA", 4) | 
                is_outlier(seurat.list[[i]], "percent.mt", 4) |
                is_outlier(seurat.list[[i]], "percent.ribo", 4)
    seurat.list[[i]]$outlier <- outlier
    # doublets
    seurat.list[[i]] <- runscDblFinder(seurat.list[[i]], assay='RNA')
    # save result
    saveRDS(seurat.list[[i]], paste0('../data/figure1/scRNA-seq/QC/',sample_names[i],'.rds'))
}

# remove outliers and doublets

dir='../data/figure1/scRNA-seq/QC/'

skin=loadRDSData(dir, sample_names, merge=TRUE)
saveRDS(skin, '../data/figure1/scRNA-seq/HF_RNA_HF_RNA_raw.rds')


skin = subset(x = skin, subset = nCount_RNA > 0)
outlier <- is_outlier(skin, "nCount_RNA", 4) | 
           is_outlier(skin, "nFeature_RNA", 4) | 
           is_outlier(skin, "percent.mt", 4) |
           is_outlier(skin, "percent.ribo", 4)
skin$outlier <- outlier

skin <- subset(skin, outlier==FALSE)
skin <- subset(skin, scDblFinder.class=='singlet')

# get blacklist genes
blacklist.genes=GetBlackListGenes(skin, MT=TRUE, Ribo=TRUE, Cellcycle=FALSE, Sex=FALSE)

saveRDS(skin,'../data/figure1/scRNA-seq/HF_RNA_postQC_merge.rds')

# load merge object post quality control
obj = readRDS('../data/figure1/scRNA-seq/HF_RNA_postQC_merge.rds')

# Integration

batch = unique(obj$batch)
batch1 = lapply(batch, function(x){
    x = str_split(x, pattern = '-') %>% unlist()
    x = x[1:3]
    x = paste0(x, collapse = '-')
    return(x)
}) %>% unlist()

meta = data.frame(row.names = batch, batch1 = batch1)
obj$batch1 = meta[obj$batch, 'batch1']

Idents(obj) = obj$seurat_clusters

HF.list = SplitObject(obj, split.by = 'batch1')
HF.list = lapply(HF.list, function(x){
  DefaultAssay(x) = 'RNA'
  x = x %>% NormalizeData(., verbose = F) %>% FindVariableFeatures(., verbose = F) %>% ScaleData(., verbose = F) %>% RunPCA(., npcs = 50, verbose = F)
})

features <- SelectIntegrationFeatures(object.list = HF.list)
anchors <- FindIntegrationAnchors(object.list = HF.list, anchor.features = features, reduction = "rpca")
HF.combined <- IntegrateData(anchorset = anchors)

DefaultAssay(HF.combined) <- "integrated"

HF.combined <- ScaleData(HF.combined, verbose = FALSE)
HF.combined <- RunPCA(HF.combined, npcs = 30, verbose = FALSE, reduction.name = 'inter_pca')
HF.combined <- RunUMAP(HF.combined, reduction = "inter_pca", dims = 1:30)
HF.combined <- FindNeighbors(HF.combined, reduction = "inter_pca", dims = 1:30)
HF.combined <- FindClusters(HF.combined)
saveRDS(HF.combined, file = '../data/figure1/scRNA-seq/HF_RNA_postQC_rpca.rds')

# Plot

colorlist = read.delim('../data/RCTD.color_LastlEdition240612.txt')
colors = colorlist$Color
names(colors) = colorlist$order

age.color = c('#3c1686', '#7696ca', '#99c9ec', '#5ea999', '#2e7737')

sample.color = c('#99CCFF', '#846DB1', '#F3DCBB')
names(sample.color) = c('HB', 'AB', 'AF')

anno = read.csv('../data/figure1/scRNA-seq/annotation.csv')
HF.combined$celltype_1 = anno[match(HF.combined$seurat_clusters, anno$seurat_clusters),]$celltype
p1 = DimPlot(HF.combined, group.by = 'celltype_1', pt.size=0.01, raster = F, cols = colors, shuffle = T)+NoLegend()+
    theme_void() + 
    theme(panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          plot.title = element_blank(),
          legend.background = element_rect(fill = "black"),
          legend.title = element_text(color = "white"),
          legend.text = element_text(color = "white"))

unique(HF.combined$batch1)
age = HF.combined$batch1 %>% str_split(., pattern = '-')# %>% unlist()
age = lapply(age, function(age){
    age = age[3]
    age = gsub('y', '', age) %>% as.numeric()
    return(age)
}) %>% unlist()

HF.combined$age_group = case_when(HF.combined$age >=3 & HF.combined$age <13 ~ "Child",
                                  HF.combined$age >=13 & HF.combined$age <18 ~ "Teenager",
                                  HF.combined$age >=18 & HF.combined$age <40 ~ "Young",
                                  HF.combined$age >=40 & HF.combined$age <60 ~ "Mid-age",
                                  HF.combined$age >=60 ~ "Old")
HF.combined$age_group = factor(HF.combined$age_group, ordered = T, levels = c('Child', 'Teenager', 'Young', 'Mid-age', 'Old'))

# sample_type
HF.combined$sample_type = ""
HF.combined@meta.data[grep('AF', HF.combined$batch), ]$sample_type = 'AF'
HF.combined@meta.data[grep('AB', HF.combined$batch), ]$sample_type = 'AO'
HF.combined@meta.data[grep('HB', HF.combined$batch), ]$sample_type = 'HO'
HF.combined$sample_type = factor(HF.combined$sample_type, ordered = T, levels = c('HO','AO','AF'))

# age_group
Idents(HF.combined) = 'sample_type'
HO = subset(HF.combined, idents = 'HO')
p2 = DimPlot(HO, group.by = 'age_group', pt.size=0.01, raster = F, cols = age.color, shuffle = T)+NoLegend()+
    theme_void() + 
    theme(panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          plot.title = element_blank(),
          legend.background = element_rect(fill = "black"),
          legend.title = element_text(color = "white"),
          legend.text = element_text(color = "white"))

# sample_type
Idents(HF.combined) = 'sample_type'
AGA = subset(HF.combined, idents = c('AF','AO'))
p3 = DimPlot(AGA, group.by = 'sample_type', pt.size=0.01, raster = F, cols = sample.color, shuffle = T)+NoLegend()+
    theme_void() + 
    theme(panel.background = element_rect(fill = "black"),
          plot.background = element_rect(fill = "black"),
          plot.title = element_blank(),
          legend.background = element_rect(fill = "black"),
          legend.title = element_text(color = "white"),
          legend.text = element_text(color = "white"))

sampleid = read.csv('../data/scRNA_Spatial_sample_NewName.csv',row.names=1)

HF.combined$batch2 = sampleid[match(HF.combined$batch, sampleid$scRNA),]$scRNA_new

Cellnum <- table(HF.combined$batch2)
Cellnum <- as.data.frame(Cellnum)

colnames(Cellnum) <- c("sample","Freq")
Cellnum$sample_type <- gsub("\\d+", "", sapply(strsplit(as.character(Cellnum$sample), "\\."), "[[", 1)) %>% factor(., levels = c('HO','AO','AF'))
Cellnum$age <- gsub("y", "", sapply(strsplit(as.character(Cellnum$sample), "\\."), "[[", 2)) %>% as.numeric

Cellnum <- Cellnum %>% arrange(sample_type, age)

p4 <- ggplot(Cellnum, aes(x = sample, y = Freq, group = sample_type, fill = sample_type)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = sample.color) +
      theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 12, color = "white"), 
        axis.text.y = element_text(size = 12, color = "white"), 
        axis.line = element_blank(),
        axis.ticks.x = element_line(color = "white", size = 0.5),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "black", color = NA),
        legend.background = element_rect(fill = "black"),
        legend.text = element_text(color = "white",size = 12))

saveRDS(HF.combined, file = '../data/figure1/scRNA-seq/HF_RNA_postQC_rpca.rds')
