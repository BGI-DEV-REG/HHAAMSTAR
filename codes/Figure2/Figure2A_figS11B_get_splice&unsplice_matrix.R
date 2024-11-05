library(dplyr)
library(data.table)
library(Matrix)
library(Seurat)
library(ggplot2)
library(SeuratObject)
library(SeuratDisk)

source('../data/figure1/Spatial/02.Cluster/LoadBGI_Spatial_fish.R')
source('../data/figure1/Spatial/02.Cluster/LoadBGI_Spatial_splice_Hair_follicle.R')

wd = paste0(args$input)
rd = paste0(args$RCTD)
setwd(wd)

input_path = paste0('../data/figure1/Spatial/00.data/',id)
output_path = paste0('../data/figure2/Spatial/04.RCTD/dynamo/',sample,'/',id)
anno_path = paste0('../data/figure1/Spatial/04.RCTD/',id)
setwd(input_path)

####################导入unsplice的数据######################
obj = LoadBGI_Spatial(paste0(sample,'_scgem.csv.gz'),
                      outdir = output_path,
                      bin_data = F,
                      bin_size = 50,
                      cell_mask = T,
                      area_mask = F,
                      save_as = "rds",
                      #pro_name = paste0(id,"_unsplice"),
                      pro_name = id,
                      UMI_GreyScale_Image = F,
                      assay = "Spatial",
                      slice = id,
                      delete_bg = T,csv=F)

obj <- obj[,obj$nCount_Spatial > 50] 

obj_s = LoadBGI_Spatial_splice(paste0(sample,'_scgem.csv.gz'),
                      outdir = output_path,
                      bin_data = F,
                      bin_size = 50,
                      cell_mask = T,
                      area_mask = F,
                      save_as = "rds",
                      pro_name = paste0(id,"_splice"),
                      UMI_GreyScale_Image = F,
                      assay = "Spatial",
                      slice = id,
                      delete_bg = T,csv=F)

obj_s = obj_s[,Cells(obj_s) %in% rownames(obj@meta.data)]

# load data with annotation
anno <- readRDS(paste0(rd,"/",id,".SCT.RCTD.rds"))
celltype <- anno$RCTD
names(celltype) <- rownames(anno@meta.data)

obj = obj[, Cells(anno)]
obj_s = obj_s[, Cells(anno)]
obj$celltype <- celltype[rownames(obj@meta.data)]

unsplice = obj@assays$Spatial@counts - obj_s@assays$Spatial@counts

# Convert as h5ad file
library(reticulate)
use_python("/opt/conda/bin/python3.8")
genes <- as.data.frame(rownames(obj), row.names = rownames(obj))
names(genes) <- "Gene"

cells <- data.frame(cell = colnames(obj),celltype = obj$celltype, row.names = colnames(obj))
names(cells) <- c("CellID","celltype")

row <- obj@images[[1]]@coordinates$row
col <- obj@images[[1]]@coordinates$col
coordinates <- list(matrix(c(row, col), ncol = 2))
names(coordinates) <- "spatial"

ann <- import("anndata")
ad <- ann$AnnData(X = obj@assays$Spatial@data, obs = genes, var = cells, varm = coordinates, layers = list(counts = obj@assays$Spatial@counts,spliced = obj_s@assays$Spatial@counts, unspliced = unsplice))
ad <- ad$T

ad$write_h5ad(file.path(output_path, paste0(id, "_all_assays.h5ad")))
