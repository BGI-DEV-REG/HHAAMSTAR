rm(list=ls())
gc()

library(ggplot2)
library(dplyr)
library(stringr)
library(patchwork)

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

spatial.meta = read.csv('../data/figure1/Spatial/03.Spatial_cluster/A03090A4_9HB-16y-M_SCT_selection_b.Spatial.metadata_0807.csv', row.names = 1)
RCTD.meta = read.csv('../data/figure1/Spatial/03.Spatial_cluster/A03090A4_9HB-16y-M_SCT_selection_b.SCT.metadata.csv')
rownames(RCTD.meta) = RCTD.meta$cellid

spatial.meta$celltype = RCTD.meta[spatial.meta$CellID, 'celltype']

percentage = spatial.meta %>% select(celltype, spatial_clusters) %>%
  filter(celltype != 'unknown') %>%
  group_by(celltype, spatial_clusters) %>%
  summarise(n = n()) %>%
  mutate(percentage = n/sum(n))

new_order = c('Mitotic zone',"Elongation zone","Keratogenouse zone",'Sloughing/Consolidation zone','Bulge','Follicle neck','Epidermis')
percentage$spatial_clusters = factor(percentage$spatial_clusters, levels = new_order%>%rev)

colorlist = colorlist[c(1:length(new_order))]
names(colorlist) = sort(new_order)

ct_order = percentage %>%
    filter(spatial_clusters == 'Mitotic zone') %>%
    arrange(percentage) %>%
    pull(celltype)

options(repr.plot.width = 8, repr.plot.height = 12)
pdf('Figure1E_RCTD-SpatialCluster_pct.pdf', width = 8, height = 12)
ggplot(percentage, aes(x = celltype, y= percentage, fill = spatial_clusters,
                       stratum = spatial_clusters, alluvium = spatial_clusters)) +
  geom_bar(stat="identity")+
  theme_classic() +
  labs(x='percentage',y = 'celltype', title = 'celltype proporation in spatial_clusters')+
  coord_flip()+
  scale_fill_manual(values = colorlist)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(panel.border = element_blank())+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 10, color = "white", angle = 0,vjust = 0.5,hjust = 0.5),
        axis.text.y = element_text(size = 10, color = "white"),
        axis.line = element_line(color = "white", size = 0.5),
        axis.ticks.x = element_line(color = "white", size = 0.5),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),#element_rect(fill = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "black", color = NA),
        legend.background = element_rect(fill = "black"),
        legend.text = element_text(color = "white",size = 10))+
scale_x_discrete(limits = ct_order)
dev.off()
