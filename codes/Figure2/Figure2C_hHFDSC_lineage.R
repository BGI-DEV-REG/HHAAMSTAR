rm(list=ls())
gc()

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(grid)
library(tibble)
library(stringr)
library(parallel)
library(RColorBrewer)
library(future)

options(future.globals.maxSize = 200 * 1024^3)
setwd('../data/figure2/monocle3/HFDSC')

colors = read.delim('../data/RCTD.color_LastlEdition240612.txt')
names = colors$order
colors = colors$Color
names(colors) = names

# subset hHFDSC-lineage related celltypes and re-integrated

obj=readRDS('../data/figure1/scRNA-seq/HF_RNA_postQC_rpca.rds')
Idents(obj) = obj$celltype

HF.obj = subset(obj, idents = c('hfDSC','Dermal_papilla','Dermal_sheath'))
Idents(HF.obj) = HF.obj$batch1
HF.obj = subset(HF.obj, ident = '8-HB-45y', invert = T) #8-HB-45y 22cells
DefaultAssay(HF.obj) = 'RNA'

HF.list = SplitObject(HF.obj, split.by = 'batch1')

HF.list = lapply(HF.list, function(x){
    DefaultAssay(x) = 'RNA'
    x = x %>% NormalizeData(., verbose = F) %>% FindVariableFeatures(., verbose = F) %>% ScaleData(., verbose = F) %>% RunPCA(., npcs = 50, verbose = F)
})
features <- SelectIntegrationFeatures(object.list = HF.list)
anchors <- FindIntegrationAnchors(object.list = HF.list, anchor.features = features, reduction = "cca") # k.filter = 80
obj <- IntegrateData(anchorset = anchors, k.weight = 80)

DefaultAssay(obj) <- "integrated"
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE, reduction.name = 'inter_pca')
seed = 42
dist = 1.3
obj <- RunUMAP(obj, reduction = "inter_pca", dims = 1:30, seed.use = seed, verbose = F, min.dist = dist)

saveRDS(obj, file = '../data/figure2/monocle3/HFDSC/HF_RNA_HFDSC_lineage_cca.rds')

p = DimPlot(obj, group.by = 'celltype', cols = colors, pt.size = 0.5)+NoLegend()+
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              plot.title = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank())
p
ggsave(p, file = paste0('FigS3_dimplot.png'), width = 5, height = 5)

gene = c('RBP4','HAPLN1','SOX2','COL11A1','LEF1','CDK1')
options(repr.plot.width = 30, repr.plot.height = 5)
lapply(gene, function(x){
    p = FeaturePlot(obj, features = x, cols = viridis::viridis(5)[-2], pt.size = 0.5, raster = F, order = T)+NoLegend()+
        theme(axis.text = element_blank(),
              axis.title = element_blank(),
              plot.title = element_blank(),
              axis.line = element_blank(),
              axis.ticks = element_blank()
              )
    ggsave(p, file = paste0('FigS3_',x, '_featureplot.png'), width = 3, height = 3)
})



# monocle3 analysis

library(monocle3)

DefaultAssay(obj) = 'RNA'

data <- GetAssayData(obj, assay = 'RNA', slot = 'counts')
cell_metadata <- obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
              
## pre-process
cds <- preprocess_cds(cds, num_dim = 50)
# plot_pc_variance_explained(cds)

cds <- reduce_dimension(cds,preprocess_method = "PCA")
# plot_cells(cds)

cds <- cluster_cells(cds, reduction_method = "UMAP")

fData(cds)$gene_short_name <- rownames(fData(cds))

recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)
head(recreate.partitions)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

list.cluster <-obj@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster

cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- obj@reductions$umap@cell.embeddings
cds <- learn_graph(cds, use_partition = TRUE)

## set root cell
get_earliest_principal_node <- function(cds, time_bin=c('DEGC')){
  # 
  cell_ids <- which(colData(cds)[, "celltype"] == time_bin)
  # 
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  # 
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

nodes_vec <- get_earliest_principal_node(cds,"hfDSC")
cds = order_cells(cds, root_pr_nodes=nodes_vec,reduction_method = "UMAP")

cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

options(repr.plot.width = 10, repr.plot.height = 10)
p1 = plot_cells(cds, 
                reduction_method="UMAP", 
                color_cells_by="celltype", 
                trajectory_graph_color = "white", 
                show_trajectory_graph = T, 
                label_leaves = F, 
                label_roots = F,
                cell_size = 1,
                trajectory_graph_segment_size = 1.5,
                label_branch_points = F, 
                label_cell_groups = F) + 
    scale_color_manual(values = colors)+
    theme_void() + 
    theme(legend.position = "none",
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank())
p1
ggsave(plot = p1, file = paste0("cds_pseudotime_Celltype.pdf"), width = 10, height = 10)
ggsave(plot = p1, file = paste0("cds_pseudotime_Celltype.png"), width = 10, height = 10)

p2 = plot_cells(cds = cds,
                color_cells_by = "pseudotime",
                trajectory_graph_color = "white", 
                show_trajectory_graph = T, 
                label_leaves = F, 
                label_roots = F,
                cell_size = 1,
                trajectory_graph_segment_size = 1.5,
                label_branch_points = F, 
                label_cell_groups = F) + viridis::scale_color_viridis(option = "D") + 
    theme_void() + 
    theme(legend.position = "right",
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.line = element_blank())
p2
ggsave(plot = p2, file = paste0("cds_pseudotime.pdf"), width = 10, height = 10)
ggsave(plot = p2, file = paste0("cds_pseudotime.png"), width = 10, height = 10)

## trajectory genes enrichment
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=8)
Track_genes <- Track_genes[,c(2,3,4,5,6)] %>% filter(q_value < 1e-3)
write.csv(Track_genes, "Trajectory_genes.csv", row.names = F)

### heatmap
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

Track_genes1 = Track_genes %>% filter(p_value < 0.05) %>% top_n(3000, wt = morans_I)
dim(Track_genes1)

pdata = pData(cds)
obj$pseudotime = pdata[Cells(obj), 'monocle3_pseudotime']

rownames(abbreviation) = colorlist$order
obj$celltype_2 = colors[obj$celltype, 'abbreviation']
unique(obj$celltype_2) %>% sort

lineage = c('hHFDSC','DP','DS')
obj$celltype_2 = factor(obj$celltype_2, level = lineage)
levels(obj$celltype_2)

meta = obj@meta.data
meta = meta %>% tibble::rownames_to_column(.,'Cell') %>% 
    group_by(celltype_2) %>% 
    arrange(celltype_2,pseudotime) %>% 
    mutate(bin = ntile(pseudotime, 20)) %>% 
    mutate(rank = paste0(as.character(celltype_2),'.',bin)) %>%
    tibble::column_to_rownames(.,'Cell')
head(meta,3)

obj$bin = meta[Cells(obj), 'bin']
obj$rank = meta[Cells(obj), 'rank']
head(obj,3)

exp <- t(as.matrix(GetAssayData(obj, assay = 'RNA', slot = 'data')))
HFSC_exp = aggregate(exp, FUN = mean, by = list(obj$rank))
HFSC_exp %>% 
  column_to_rownames("Group.1") -> df

df_t <- t(df)
df_t_gene <- df_t[rownames(df_t) %in% traj_gene1$gene_short_name,]
gene_sd<-apply(df_t_gene,1,sd)
df_sd <- cbind(df_t_gene,gene_sd)
df_sd[1:3,]
dim(df_sd)
write.csv(df_sd,"DSC_sd.csv")

df_sd = read.csv('DSC_sd.csv',row.names=1)

df_sd = df_sd[, -ncol(df_sd)]
df_sd_t <- t(df_sd)
head(df_sd_t,3)

mat_cha <- t(scale(df_sd_t))
head(mat_cha,3)

head(as.data.frame(mat_cha),3)

order = c(unique(meta$rank)[grep('hHFDSC', unique(meta$rank))],
          unique(meta$rank)[grep('DP', unique(meta$rank))],
          unique(meta$rank)[grep('^DS', unique(meta$rank))])
order

mat_cha[mat_cha >= 3] = 3
mat_cha[mat_cha <= -3] = -3
mat_cha = mat_cha[, order]
mat_cha_smooth <- mat_cha[, order]
mat_cha_smooth <- t(apply(mat_cha_smooth,1,function(x){smooth.spline(x,df=3)$y}))
mat_cha_smooth[mat_cha_smooth >= 3] = 3
mat_cha_smooth[mat_cha_smooth <= -3] = -3
head(mat_cha_smooth,3)

colnames(mat_cha_smooth) = order

set.seed(2024)

options(repr.plot.width = 5, repr.plot.height = 5)
k_means <- 15
p_clust=pheatmap::pheatmap(mat_cha,
                     cutree_rows=k_means, breaks = seq(-3,3, by = 0.1),
                     cluster_cols = F,show_rownames = FALSE,
                     show_colnames = FALSE)
split_matrix=data.frame(cutree(p_clust$tree_row,k=k_means))

data.meta = meta %>% select(celltype_2, rank) %>% distinct
data.meta$color = colors[match(data.meta$celltype_2, colors$abbreviation), 'Color']

column_ha = HeatmapAnnotation(
  pattern = data.meta$rank,
  col = list(pattern = setNames(data.meta$color, data.meta$rank))
)

colnames(mat_cha_smooth)

options(repr.plot.width = 15, repr.plot.height = 18)
clust_hm2 <- Heatmap(mat_cha_smooth,name='Exp',
                     circlize::colorRamp2(c(-3,0,0.2,1,2), brewer.pal(11, "RdBu")[c(10,6,4,2,1)]),
                     # colorRamp2(c(-3,0,3), c('#3E9BFEFF','white','#F05B12FF')),
                     # colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(100),
                     left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:21),
                                                                      labels_gp = gpar(col = "white", fontsize = 8))),
                     row_split = split_matrix,
                     top_annotation = column_ha,
                     cluster_columns = FALSE,
                     #right_annotation = ha,
                     show_row_names = FALSE,
                     show_column_names = FALSE)
clust_hm2

split_matrix1 = split_matrix %>% 
    mutate(order = p_clust$tree_row$order) %>%
    # rownames_to_column(.,'gene') %>%
    rename('cutree.p_clust.tree_row..k...k_means.' = 'module') %>%
    arrange(module)
unique(split_matrix1$module)
head(split_matrix1)

module.levels = c(
    1,6,7,11,10,3,4,9,14,5,13,12,8,15,2
)%>% rev
genelist = split_matrix1
genelist$module = factor(genelist$module, level = module.levels)
genelist = genelist %>% rownames_to_column(.,'gene') %>% group_by(module) %>% arrange(module,order)
head(genelist,3)

tf = read.delim('../data/database/allTFs_hg38.txt',header=F) %>% dplyr::rename('TF' = 'V1')
head(tf,3)

scenic = read.csv('/data/work/01.HF_RNA/HFDSC_lineage/SCENIC/ctx.csv')
colnames(scenic) = c(scenic[2,c(1:2)], scenic[1,-c(1:2)])
scenic = scenic[-c(1:2),]
# head(scenic,3)

gene = c('RBP4','HAPLN1','ASPN','SOX2','LEF1','MMP11','COL11A1','COL3A1',
         'NDUFS7','UQCRH','COX6C,HK1','LDHA,CDH1','GRHL3','CLDN1,DNAJC4','YOD1',
         'EIF2AK3','ELB','FOSL2','MAP2K3','IL6R','TNFRSF1B','CXCR4','ITGA1','ITGA8','NPNT','PTCH1','GLI2','GLI3',
         'CCND1','CDK2','TOP2A','SPP1','IGF1','IGF2',
         'FAS','TNFRSF12A','GADD45A','TRAF1','TNFRSF11B','MAP2K4','APCDD1','PRICKLE1','TBX18','TP63',
         'KRT14','KRT10','CCNB1','CDC25B','CDK6','AREG','HBEGF','ERBB4') %>% unique
gene = union(gene, scenic$TF) %>% sort
gene

rownames(mat_cha_smooth1)[which(rownames(mat_cha_smooth1) %in% gene)]

which(rownames(mat_cha_smooth1) %in% gene)

options(repr.plot.width = 10, repr.plot.height = 12)
mat_cha_smooth1 = mat_cha_smooth[genelist$gene, order]
clust_hm2 <- Heatmap(mat_cha_smooth1, name='Exp',
                     # row_order = rownames(genelist),
                     circlize::colorRamp2(c(-3,0,0.2,0.6,1.5), brewer.pal(11, "RdBu")[c(10,6,4,2,1)]),
                     # left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = 2:19),
                     #                                                  labels_gp = gpar(col = "white", fontsize = 8))),
                     # colorRamp2(c(-3,0,3), c('#3E9BFEFF','white','#F05B12FF')),
                     # left_annotation = row_ha,
                     # cluster_rows = FALSE,
                     cluster_row_slices = TRUE,
                     clustering_distance_rows = "euclidean",
                     clustering_method_rows = "complete",
                     row_split = genelist$module,
                     top_annotation = column_ha,
                     cluster_columns = FALSE,
                     # cluster_rows = FALSE,
                     #right_annotation = ha,
                     show_row_names = FALSE,
                     show_column_names = FALSE,
                     use_raster = T)+
  rowAnnotation(link = anno_mark(at = which(rownames(mat_cha_smooth1) %in% gene), 
                                   labels = rownames(mat_cha_smooth1)[which(rownames(mat_cha_smooth1) %in% gene)], labels_gp = gpar(fontsize = 6)))
clust_hm2

pdf('HF_DSC_lineage_heatmap.pdf', width = 10, height = 12)
clust_hm2
dev.off()