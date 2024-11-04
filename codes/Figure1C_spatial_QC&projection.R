library(dplyr)
library(data.table)
library(Matrix)
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(SeuratDisk)

source('../data/figure1/code/LoadBGI_Spatial.R')

setwd(input_path)
input_path = paste0('../data/figure1/Spatial/raw_data/',id)
nCount_Spatial = 50
output_path = paste0('../data/figure1/Spatial/QC/',id)
result_path = paste0('../data/figure1/Spatial/RCTD/',id)

#### Creat obj 
obj = LoadBGI_Spatial(paste0(id,'_scgem.csv.gz'),
                      outdir = getwd(),
                      bin_data = F,
                      bin_size = 50,
                      cell_mask = T,
                      area_mask = F,
                      save_as = "rds",
                      pro_name = id,
                      UMI_GreyScale_Image = F,
                      assay = "Spatial",
                      slice = id,
                      delete_bg = T,csv=F)

#### QC and selecting cells for further analysis
obj <- obj[,obj$nCount_Spatial > nCount_Spatial] # nCount_RNA can be adjusted

#### SCT
DefaultAssay(obj) = 'Spatial'
obj <- SCTransform(obj, verbose = FALSE,assay = 'Spatial')
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE) ## or assay = "RNA"
obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
obj <- FindClusters(obj, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30)

saveRDS(obj,paste0(output_path,'/',id,'.SCT.rds'))

#### Meta.data 
meta = as.data.frame(obj@meta.data)
write.csv(meta,paste0(output_path,'/',id,'.Seurat.metadata.SCT.csv'),quote=F)

#### Convert as h5ad file
library(reticulate)
use_python("/opt/conda/bin/python3.8")

genes <- as.data.frame(rownames(obj), row.names = rownames(obj))
names(genes) <- "Gene"

cells <- as.data.frame(colnames(obj), row.names = colnames(obj))
names(cells) <- "CellID"

row <- obj@images[[1]]@coordinates$row
col <- obj@images[[1]]@coordinates$col
coordinates <- list(matrix(c(row, col), ncol = 2))
names(coordinates) <- "spatial"

ann <- import("anndata")
ad <- ann$AnnData(X = obj@assays$SCT@data, obs = genes, var = cells, varm = coordinates,layers = list(counts = obj@assays$SCT@counts))

ad <- ad$T
ad$write_h5ad(file.path(output_path, paste0(id, "_SCT.h5ad")))


# RCTD
##传参（将时空组数据传进来）

library(data.table)
library(spacexr)
library(Seurat)
library(readr)
library(ggplot2)
library(dplyr)

colorlist = c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941",
             "#006FA6", "#A30059", "#FFE4E1", "#0000A6", "#63FFAC",
             "#B79762", "#004D43", "#8FB0FF", "#997D87", "#5A0007",
             "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#FF2F80",
             "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9",
             "#B903AA", "#DDEFFF", "#7B4F4B", "#A1C299", "#0AA6D8",
             "#00A087FF", "#4DBBD5FF", "#E64B35FF", "#3C5488FF", "#F38400",
             "#A1CAF1", "#C2B280", "#848482", "#E68FAC", "#0067A5",
             "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300",
             "#882D17", "#8DB600", "#654522", "#E25822", "#2B3D26",
             "#191970", "#000080",
             "#6495ED", "#1E90FF", "#00BFFF", "#00FFFF", "#FF1493",
             "#FF00FF", "#A020F0", "#63B8FF", "#008B8B", "#54FF9F",
             "#00FF00", "#76EE00", "#FFF68F")

obj <- readRDS(paste0(output_path, '/',id, '.SCT.rds'))
print("obj read successfully！")

## import scRNA-seq data
sc=readRDS("../data/figure1/scRNA-seq/HF_RNA_postQC_rpca.rds")
print("sc read successfully！")

Idents(sc) = sc$celltype
sc=subset(sc, downsample = 1000)

## get spatial location
exp_spatial=as.matrix(obj@assays$Spatial@counts)
exp_spatial=as.data.frame(exp_spatial)

obj$coor_y = obj@images[[names(obj@images)]]@coordinates$row
obj$coor_x = obj@images[[names(obj@images)]]@coordinates$col
coord <- data.frame(stereo_1 = obj$coor_x, stereo_2 = obj$coor_y)
obj[["stereo"]] <- SeuratObject::CreateDimReducObject(embeddings = as.matrix(coord), key = "stereo_", assay = "SCT")

row <- obj$coor_y
col <- obj$coor_x
coord_spatial <- data.frame(row,col)
rownames(coord_spatial)<- Cells(obj)
nUMI_spatial=obj@meta.data[,"nCount_Spatial"]
names(nUMI_spatial)=rownames(obj@meta.data)

## get scRNA-seq data
exp_sc=as.matrix(sc@assays$RNA@counts)
exp_sc=as.data.frame(exp_sc)

celltype_sc=sc$celltype
names(celltype_sc)=rownames(sc@meta.data)
celltype_sc=as.factor(celltype_sc)
nUMI_sc=sc@meta.data[,"nCount_RNA"]
names(nUMI_sc)=rownames(sc@meta.data)

## make query and reference for RCTD
reference <- Reference(exp_sc, celltype_sc, nUMI_sc)
puck <- SpatialRNA(coord_spatial, exp_spatial, nUMI_spatial)
print("data have praperad successfully！")

rm(sc)
gc()

## run RCTD
myRCTD <- create.RCTD(puck, reference, max_cores = 1,UMI_min_sigma = args$UMI)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')


setwd(result_path)

rds_name<-paste0(id, "_myRCTD.rds")
saveRDS(myRCTD, file = rds_name)
print("RCTD successfully！")

# result
head(myRCTD@results$results_df)
df = myRCTD@results$results_df
df = df[df$spot_class != 'reject', ]
table(df$first_type)

# add meta
obj$RCTD = ''
for (i in rownames(df)) {
  obj@meta.data[i, ]$RCTD = as.character(df[i, ]$first_type)
}
obj$RCTD[obj$RCTD == ""] = "unknown"

saveRDS(obj, file = paste0(id, '.SCT.RCTD.rds'))
write_csv(obj@meta.data, file = paste0(id, '.SCT.metadata.csv'))

