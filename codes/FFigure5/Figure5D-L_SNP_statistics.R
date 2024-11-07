library(MotifDb)

library(Seurat)
library(Signac)
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(rtracklayer)
library(parallel)
library(viridis)
library(patchwork)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(ComplexHeatmap)
sessionInfo()

setwd('../data/figure5/ldsc_SNP/')

# lead SNP
files = list.files('../data/figure5/ldsc_SNP/GWAS_leadsnp/',full.names=T)

leadsnp = lapply(files, function(x){
    df = read.delim(x)
}) %>% Reduce(rbind, .)
head(leadsnp,3)

leadsnp = leadsnp %>% select(uniqID, rsID)
colnames(leadsnp) = c("hg19_uniqID","rsID")
head(leadsnp,3)

write.table(leadsnp,file = 'leadsnp17.txt',quote=F,sep = '\t',row.names=F)

# liftover : hg19 - hg38
ref = read.delim('../MACSzRm2XoZXhhR2.txt') %>% select(X.Uploaded_variation,Location,REF_ALLELE,Allele) %>% distinct
colnames(ref)[1] = 'rsID'
head(ref,3)

leadsnp$Chr = sapply(ref[match(leadsnp$rsID, ref$rsID),]$Location, function(x){
    x = unlist(strsplit(x,':'))[1] %>% paste0('chr', .)
}) %>% unlist
leadsnp$Start = sapply(ref[match(leadsnp$rsID, ref$rsID),]$Location, function(x){
    x = unlist(strsplit(x,'-'))[1] %>% gsub('\\d+:','',.) %>% as.numeric
}) %>% unlist
leadsnp$hg38_REF = ref[match(leadsnp$rsID, ref$rsID),]$REF_ALLELE
leadsnp$hg38_ALT = ref[match(leadsnp$rsID, ref$rsID),]$Allele
head(leadsnp,3)

leadsnp$hg19_uniqID = gsub(':','_',leadsnp$hg19_uniqID) %>% paste0('chr', .)
head(leadsnp,3)

length(leadsnp$rsID %>% unique)
# 2021


## filter peaks (169/2021)

leadsnp.gr = GRanges(seqnames = leadsnp$Chr,
                 ranges = IRanges(leadsnp$Start, leadsnp$Start))

peaks = readRDS('../peaks_grange.rds') # extract from scATAC-seq
overlaps = findOverlaps(leadsnp.gr, peaks, type = "within")

leadsnp.peaks = leadsnp[queryHits(overlaps),]
dim(leadsnp.peaks)
leadsnp.gr.peaks = leadsnp.gr[queryHits(overlaps),]
head(leadsnp.gr.peaks,3)

peaks.leadsnp = peaks[subjectHits(overlaps),]
peaks.leadsnp

peaks.df = as.data.frame(peaks.leadsnp)
peaks.df$peak = paste0(peaks.df$seqnames,'-',peaks.df$start,'-',peaks.df$end)
head(peaks.df,3)

dim(peaks.df)

leadsnp.peaks$peak = peaks.df$peak
head(leadsnp.peaks,3)

leadsnp.peaks1 = distinct(leadsnp.peaks)
dim(leadsnp.peaks1)


## filter DEP (67/169)

dep = read.csv('../AGA_health_DEseq2_240717.csv', row.names = 1) %>% filter(abs(log2FoldChange) > 0.25 & pvalue < 0.05)
colnames(dep)[1] = 'peak'
head(dep,3)

dep1 = dep %>% select(peak) #%>% distinct
dep1$chr = sapply(dep1$peak, function(x){ unlist(strsplit(x,'-'))[1] }) %>% unlist
dep1$start = sapply(dep1$peak, function(x){ unlist(strsplit(x,'-'))[2] %>% as.numeric }) %>% unlist
dep1$end = sapply(dep1$peak, function(x){ unlist(strsplit(x,'-'))[3] %>% as.numeric }) %>% unlist
# head(dep1,3)
dep.gr = GRanges(seqnames = dep1$chr,
                 ranges = IRanges(dep1$start, dep1$end))
head(dep.gr,3)

dep.overlaps = findOverlaps(leadsnp.gr.peaks, dep.gr, type = "within")

leadsnp.dep = leadsnp.peaks[queryHits(dep.overlaps),]
dim(leadsnp.dep)
leadsnp.gr.dep = leadsnp.gr.peaks[queryHits(dep.overlaps),]
head(leadsnp.gr.dep,3)

leadsnp.dep$peak = dep[subjectHits(dep.overlaps), ]$peak
leadsnp.dep$ATAC_log2FoldChange = dep[subjectHits(dep.overlaps), ]$log2FoldChange
leadsnp.dep$ATAC_pvalue = dep[subjectHits(dep.overlaps), ]$pvalue
leadsnp.dep$ATAC_Tag = dep[subjectHits(dep.overlaps), ]$cluster
leadsnp.dep = distinct(leadsnp.dep)
head(leadsnp.dep,3)

unique(leadsnp.dep$rsID) %>% length

leadsnp.dep[which(is.na(leadsnp.dep$rsID)), ]


## p2g (55/67)

links = read.csv('../HF_ATAC_p2g.links_500kb.csv')
head(links,3)

links = links %>% select(peak, gene) %>% distinct()
head(links,3)

leadsnp.p2g = left_join(leadsnp.dep, links, by = 'peak')
dim(leadsnp.p2g)
head(leadsnp.p2g,3)

unique(leadsnp.p2g$rsID) %>% length

unique(leadsnp.p2g[leadsnp.p2g$gene != 'NA', ]$rsID) %>% length

leadsnp.p2g = leadsnp.p2g[leadsnp.p2g$gene != 'NA', ]



## DEG (52/55)

deg = readRDS('../DESeq2_0606.rds') %>% filter(abs(log2FoldChange) > 0.25 & pvalue < 0.05)
head(deg,3)

deg$RNA_Tag = paste0(deg$ident.1,'_',deg$ident.2,'.',deg$celltype)
deg = deg %>% select(gene, log2FoldChange, pvalue, RNA_Tag) %>% distinct
colnames(deg) = c('gene', 'RNA_log2FoldChange', 'RNA_pvalue', 'RNA_Tag')
head(deg,3)

leadsnp.p2g.deg = left_join(leadsnp.p2g, deg, by = 'gene') %>% distinct
dim(leadsnp.p2g.deg)
head(leadsnp.p2g.deg,3)

leadsnp.p2g.deg = leadsnp.p2g.deg[!which(is.na(leadsnp.p2g.deg$rsID)), ]

unique(leadsnp.p2g.deg$rsID) %>% length
unique(leadsnp.p2g.deg$rsID)

write.table(leadsnp.p2g.deg, file = 'Summary_snp_dp_2g_55.txt',quote = F,sep = '\t', row.names = F)

leadsnp.p2g.deg1 = leadsnp.p2g.deg[complete.cases(leadsnp.p2g.deg), ] %>% distinct
unique(leadsnp.p2g.deg1$rsID) %>% length
write.table(leadsnp.p2g.deg1 %>% dplyr::select(hg19_uniqID) %>% distinct, file = 'Summary_snp_dp_2g_deg_55P2G_forDeltaSVM.txt',quote = F,sep = '\t', col.names = F, row.names = F)

## delta SVM

leadsnp.p2g.deg = read.delim('Summary_snp_dp_2g_55.txt')
head(leadsnp.p2g.deg,3)

leadsnp.p2g.deg = leadsnp.p2g.deg %>% rename(peak = 'peakID')
head(leadsnp.p2g.deg,3)

tf = read.csv('snpsummary.pred_tf.csv',row.names=1) %>% select(snp,tf,preferred_allele) %>% rename('snp' = 'hg19_uniqID')
head(tf,3)

tf[tf$hg19_uniqID %in% leadsnp.dep$hg19_uniqID, ]$hg19_uniqID %>% unique %>% length

tf[tf$hg19_uniqID %in% leadsnp.p2g.deg$hg19_uniqID, ]$hg19_uniqID %>% unique %>% length

leadsnp.p2g.deg$tf = tf[match(leadsnp.p2g.deg$hg19_uniqID, tf$hg19_uniqID), 'tf']
head(leadsnp.p2g.deg,3)

leadsnp.p2g.deg1 = leadsnp.p2g.deg[!is.na(leadsnp.p2g.deg$tf), ]
head(leadsnp.p2g.deg1,3)


# SNP statistics 

## 
#2021-169

SNPinPeaks = data.frame(Var1 = c('beyond_peaks','in_peaks'),
                        Freq = c(1852, 169),
                        Pct = c(1852/2021 * 100, 169/2021 * 100))
SNPinPeaks

options(repr.plot.width = 6, repr.plot.height = 4)
p1 = ggplot(SNPinPeaks, aes(x = "", y = Pct, fill = Var1)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    labs(fill = "region") +
    theme_void()+
    scale_fill_manual(values = GSE_color[c(3,2)]) + 
    geom_text(aes(#label = paste0(round(Pct,2), "%")), 
              label = Freq),
              position = position_stack(vjust = 0.5), 
              color = "white") + # 设置标注颜色
    theme(plot.background = element_rect(fill = "black"),
          panel.background = element_blank(),# element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"),
          legend.background = element_blank(),
          plot.title = element_text(color = "white"),
          legend.key = element_blank())
p1
ggsave(p1, file = 'LeadSNP_in_peaks.pdf', width = 6, height = 5)

#169-67

leadSNPinDEP = data.frame(Var1 = c('SNPs_in_peaks','SNPs_in_DEPs'), 
                          Freq = c(102,67), 
                          Pct = c(102/169 * 100, 67/169 * 100))
leadSNPinDEP

options(repr.plot.width = 6, repr.plot.height = 4)
p2 = ggplot(leadSNPinDEP, aes(x = "", y = Pct, fill = Var1)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    labs(fill = "region") +
    theme_void()+
    scale_fill_manual(values = GSE_color[c(1,2)]) + 
    geom_text(aes(#label = paste0(round(Pct,2), "%")), 
              label = Freq),
              position = position_stack(vjust = 0.5), 
              color = "white") + # 设置标注颜色
    theme(plot.background = element_rect(fill = "black"),
          panel.background = element_blank(),# element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"),
          legend.background = element_blank(),
          plot.title = element_text(color = "white"),
          legend.key = element_blank())
p2
ggsave(p2, file = 'LeadSNP_in_DEPs.pdf', width = 6, height = 5)

#67-55

leadSNPlinkGene = data.frame(Var1 = c('Linked genes','No linked genes'), 
                          Freq = c(55,12), 
                          Pct = c(55/67 * 100, 12/67 * 100))
leadSNPlinkGene


options(repr.plot.width = 6, repr.plot.height = 4)
p3 = ggplot(leadSNPlinkGene, aes(x = "", y = Pct, fill = Var1)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    labs(fill = "region") +
    theme_void()+
    scale_fill_manual(values = c(GSE_color[1],'#206fce')) + 
    geom_text(aes(#label = paste0(round(Pct,2), "%")), 
              label = Freq),
              position = position_stack(vjust = 0.5), 
              color = "white") + # 设置标注颜色
    theme(plot.background = element_rect(fill = "black"),
          panel.background = element_blank(),# element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"),
          legend.background = element_blank(),
          plot.title = element_text(color = "white"),
          legend.key = element_blank())
p3
ggsave(p3, file = 'LeadSNP_in_DEPs_linked_genes.pdf', width = 6, height = 5)

#67-21

leadSNPreTF = data.frame(Var1 = c('NO TF','TF'), 
                          Freq = c(46,21), 
                          Pct = c(46/67 * 100, 21/67 * 100))
leadSNPreTF

options(repr.plot.width = 6, repr.plot.height = 4)
p3 = ggplot(leadSNPreTF, aes(x = "", y = Pct, fill = Var1)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    labs(fill = "region") +
    theme_void()+
    scale_fill_manual(values = c('#2171a8','#ee7b20')) + 
    geom_text(aes(#label = paste0(round(Pct,2), "%")), 
              label = Freq),
              position = position_stack(vjust = 0.5), 
              color = "white") + # 设置标注颜色
    theme(plot.background = element_rect(fill = "black"),
          panel.background = element_blank(),# element_rect(fill = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"),
          legend.background = element_blank(),
          plot.title = element_text(color = "white"),
          legend.key = element_blank())
p3
ggsave(p3, file = 'LeadSNP_related_TF.pdf', width = 6, height = 5)

## snp region annotation

library(clusterProfiler)
library("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak = leadsnp.dep[, c('Chr','Start','Start')] %>% distinct
colnames(peak) = c('chr','start','end')
# peak$end = peak$start+1
head(peak,3)
dim(peak)

snp67.gr = GRanges(seqnames = peak$chr,
                    ranges = IRanges(peak$start, peak$end))

library(ChIPseeker)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- getTagMatrix(snp67.gr, windows=promoter)

peakAnno <- annotatePeak(
    snp67.gr,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db")

peakAnno1 = peakAnno@annoStat %>% arrange(desc(Frequency))
peakAnno1

options(repr.plot.width = 3.5, repr.plot.height = 4)
p3 = ggplot(peakAnno1, aes(x = Frequency, y = Feature, fill = Feature)) +
  geom_bar(stat = "identity", , position = "stack") +
  scale_fill_manual(values = archrcolor[5:13]) +
  theme_minimal() +
  scale_y_discrete(limits=peakAnno1$Feature %>% rev)+
  theme(legend.position = 'none',
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "black"),
        panel.border = element_rect(fill = NA, color = "white",
                                    size = 1, linetype = "solid"),
        plot.background = element_rect(fill = "black"),
        axis.text.x = element_text(color = "white", size = 10, angle = 0, vjust=0.5, hjust=1),
        axis.text.y = element_text(color = "white", size = 10),
        legend.text = element_text(color = "white"),
        plot.title = element_text(color = "white"))
p3
ggsave(p3, file = 'SNP67_region_anno.pdf',width = 3.5, height = 4)

## SNP & TF

head(tf,3)

leadsnp.dep = left_join(leadsnp.dep,tf,by = 'hg19_uniqID') %>% distinct
head(leadsnp.dep,3)

leadsnp.dep$tf = tf[match(leadsnp.dep$hg19_uniqID, tf$hg19_uniqID), 'tf']
leadsnp.dep$preferred_allele = tf[match(leadsnp.dep$hg19_uniqID, tf$hg19_uniqID), 'preferred_allele']

snp.tf = leadsnp.dep %>% select(rsID,tf,preferred_allele) %>% distinct
snp.tf = snp.tf[complete.cases(snp.tf), ]
head(snp.tf,3)


snp.tf.sum = snp.tf %>% group_by(rsID,preferred_allele) %>% summarise(n_tf = n())
snp.tf.sum[snp.tf.sum$preferred_allele == 'Loss', ]$n_tf = -snp.tf.sum[snp.tf.sum$preferred_allele == 'Loss', ]$n_tf
snp.tf.sum

snp.tf.order = c('rs9388486','rs2093283','rs2612560','rs3778607','rs11029994','rs140371629','rs16892288','rs2298274','rs4894405','rs7621843','rs7642052','rs847162',
                 'rs199522','rs2163799','rs4784464','rs2024021','rs2237257','rs3012384','rs35624710','rs4946689','rs59189180')

options(repr.plot.width = 2.75, repr.plot.height = 3)
p = ggplot(snp.tf.sum, aes(x = n_tf, y = rsID, fill = preferred_allele)) +
  geom_bar(position = "stack", stat = "identity") +
  # theme_classic()+
  scale_fill_manual(values = c("Gain" = '#ff6347', "Loss" = '#00ffff')) +
  scale_y_discrete(limits = rev(snp.tf.order))+
  theme_minimal() + #coord_flip()+
  labs(title = 'SNP Gain and Loss Counts')+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(c(3, 2), "cm"),
        panel.background = element_rect(fill = "black"),
        panel.border = element_rect(fill = NA, color = "white",size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "black"),
        axis.text.x = element_text(color = "white", size = 6, angle = 0, vjust=0.5, hjust=0.5),
        axis.text.y = element_text(color = "white", size = 6),
        legend.text = element_text(color = "white", size = 6),
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(color = "white", size = 8, vjust=0, hjust=0.5)
  )
p
ggsave(p, file = 'SNP_TF_count.pdf',width=2.75,height=3)

tf.snp.sum = snp.tf %>% group_by(tf,preferred_allele) %>% summarise(n_snp = n())
tf.snp.sum[tf.snp.sum$preferred_allele == 'Loss', ]$n_snp = -tf.snp.sum[tf.snp.sum$preferred_allele == 'Loss', ]$n_snp
tf.snp.sum

tf.snp.order = c('UBP1','NFIC','NFIX','ATF2','ATF3','ATOH1','CREB5','HSF2','NEUROG1','SPIB','T','TFAP2E',
                 'LHX9','VAX2','CEBPE','DLX4','EMX1','HLF','HOXD10','HSF4','ISX','LHX6','LHX8','MEIS3','OVOL1','RAX','SKOR2')

options(repr.plot.width = 2.75, repr.plot.height = 3)
p = ggplot(tf.snp.sum, aes(x = n_snp, y = tf, fill = preferred_allele)) +
  geom_bar(position = "stack", stat = "identity") +
  # theme_classic()+
  scale_fill_manual(values = c("Gain" = '#ff6347', "Loss" = '#00ffff')) +
  scale_y_discrete(limits = rev(tf.snp.order))+
  theme_minimal() + #coord_flip()+
  labs(title = 'TF Gain and Loss Counts')+
  theme(panel.grid = element_blank(),
        panel.spacing = unit(c(3, 2), "cm"),
        panel.background = element_rect(fill = "black"),
        panel.border = element_rect(fill = NA, color = "white",size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "black"),
        axis.text.x = element_text(color = "white", size = 6, angle = 0, vjust=0.5, hjust=0.5),
        axis.text.y = element_text(color = "white", size = 6),
        legend.text = element_text(color = "white", size = 6),
        legend.title = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        plot.title = element_text(color = "white", size = 8, vjust=0, hjust=0.5)
  )
p
ggsave(p, file = 'TF_SNP_count.pdf',width=2.75,height=3)



# peak links gene

atac.all = readRDS('../HF_ATAC_peak2gene500kb.rds')

colors1 = read.delim('../data/RCTD.color_LastlEdition240612.txt')
rownames(colors1) = colors1$order

atac = atac.all
Idents(atac) = 'age_group'
atac = subset(atac, idents = c('Mid_age','Young'))
atac@assays$ATAC@links = atac.all@assays$ATAC@links

atac$celltype = factor(atac$celltype, levels = unique(atac$celltype) %>% sort)
levels(atac$celltype)

atac$sample = factor(atac$sample, levels = c('HO','AO','AF'))
levels(atac$sample)

snp_data = read.delim('Summary_snp_dp_2g_55.txt',header=T)
snp_data1 = snp_data %>% select(rsID,Chr,Start,peak) %>% distinct
head(snp_data1,3)

snp_select = data.frame(snp = c('rs2237257','rs847162','rs2612560'),
                        celltype = c('DP',
                                     'ORSS,MELA,TAC1,hHFDSC,DS,DP',
                                     'IRS,IRSC,ORSB,ORSS'),
                        gene = c('TRAF3IP2','TFAP2E,HOXD','SLC14A2'))
# snp_select
snp_select$peak = snp_data1[match(snp_select$snp, snp_data1$rsID), ]$peak
snp_select$Chr = snp_data1[match(snp_select$snp, snp_data1$rsID), ]$Chr
snp_select$Start = snp_data1[match(snp_select$snp, snp_data1$rsID), ]$Start
snp_select$Start1 = snp_select$Start-10
snp_select$End = snp_select$Start+10
snp_select$Location = paste0(snp_select$Chr,'-',snp_select$Start1,'-',snp_select$End)
snp_select

atac$tissue = 'hair_collicle'

n = 1
snp = snp_select$snp[n]
snp_peak = snp_select$Location[n]
celltype = snp_select$celltype[n] %>% strsplit(.,',') %>% unlist
peak = snp_select$peak[n]
location = snp_select$Location[n]
Idents(atac) <- "celltype" 

show_snp <- StringToGRanges(peak) 
color <- 'white'
mcols(show_snp)$color <- color 
Idents(atac) <- "tissue" 
plot1 <- CoveragePlot(object = atac, region = peak, extend.upstream = 100000, extend.downstream = 100000, 
                      peaks = F, links = T, annotation = T, heights = 4, ncol = 1,region.highlight = show_snp)
pdf('rs2237257_peak100k.pdf', width = 8, height = 5)
plot1 & scale_fill_manual(values = cell_color[9]) & theme_minimal()+#ggtitle(paste0(snp, '-', ct))+
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "black"),
          panel.border = element_rect(fill = NA, color = "white",size = 0.5, linetype = "solid"),
          plot.background = element_rect(fill = "black"),
          axis.line = element_line(color = "white"),
          axis.line.x = element_line(color = "white", size = 0.5),
          axis.text.x = element_text(color = "white", size = 6, angle = 0, vjust=0.5, hjust=0), 
          axis.text.y = element_text(color = "white", size = 6),
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"),
          legend.position = "right",
          plot.title = element_text(color = "white"),
          strip.text = element_text(color = "white"),
          strip.background = element_blank())
dev.off()

n = 2
snp = snp_select$snp[n]
snp_peak = snp_select$Location[n]
celltype = snp_select$celltype[n] %>% strsplit(.,',') %>% unlist
peak = snp_select$peak[n]
location = snp_select$Location[n]
Idents(atac) <- "celltype" 
show_snp <- StringToGRanges(peak) 
color <- 'white'
mcols(show_snp)$color <- color 
Idents(atac) <- "tissue" 
plot1 <- CoveragePlot(object = atac, region = peak, extend.upstream = 50000, extend.downstream = 50000, 
                      peaks = F, links = T, annotation = T, heights = 4, ncol = 1,region.highlight = show_snp)
pdf('rs847162_peak50k.pdf', width = 8, height = 4)
plot1 & scale_fill_manual(values = cell_color[9]) & theme_minimal()+  #ggtitle(paste0(snp, '-', ct))+
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "black"),
        panel.border = element_rect(fill = NA, color = "white",size = 0.5, linetype = "solid"),
        plot.background = element_rect(fill = "black"),
        axis.line = element_line(color = "white"),
        axis.line.x = element_line(color = "white", size = 0.5),
        axis.text.x = element_text(color = "white", size = 6, angle = 0, vjust=0.5, hjust=0), 
        axis.text.y = element_text(color = "white", size = 6),  
        legend.text = element_text(color = "white"),
        legend.title = element_text(color = "white"),
        legend.position = "right",
        plot.title = element_text(color = "white"),
        strip.text = element_text(color = "white"),
        strip.background = element_blank()
       )
dev.off()


n = 3
snp = snp_select$snp[n]
snp_peak = snp_select$Location[n]
celltype = snp_select$celltype[n] %>% strsplit(.,',') %>% unlist
peak = snp_select$peak[n]
location = snp_select$Location[n]
Idents(atac) <- "celltype" 
show_snp <- StringToGRanges(peak) 
color <- 'white'
mcols(show_snp)$color <- color 
Idents(atac) <- "tissue" 
plot1 <- CoveragePlot(object = atac, region = peak, extend.upstream = 500000, extend.downstream = 500000, 
                      peaks = F, links = T, annotation = T, heights = 4, ncol = 1,region.highlight = show_snp)
pdf('rs2612560_peak500k.pdf', width = 8, height = 4)
plot1 & scale_fill_manual(values = cell_color[9]) & theme_minimal()+  #ggtitle(paste0(snp, '-', ct))+
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "black"),
          panel.border = element_rect(fill = NA, color = "white",size = 0.5, linetype = "solid"),
          plot.background = element_rect(fill = "black"),
          axis.line = element_line(color = "white"),
          axis.line.x = element_line(color = "white", size = 0.5),
          axis.text.x = element_text(color = "white", size = 6, angle = 0, vjust=0.5, hjust=0),
          axis.text.y = element_text(color = "white", size = 6), 
          legend.text = element_text(color = "white"),
          legend.title = element_text(color = "white"),
          legend.position = "right",
          plot.title = element_text(color = "white"),
          strip.text = element_text(color = "white"),
          strip.background = element_blank()
         )
dev.off()


# figure5-17 snps table

snp = read.delim('Summary_snp_dp_2g_0904djx_55_1008.txt')
head(snp,3)

unique(snp$rsID) %>% length

tf = read.csv('snpsummary.pred_tf.csv',row.names=1) %>% select(snp,tf,preferred_allele) %>% rename('snp' = 'hg19_uniqID')
head(tf,3)

snp1 = left_join(snp, tf, by = c('hg19_uniqID'))
head(snp1,3)

snp2 = snp1 %>% select(rsID, peak, gene, tf, preferred_allele) %>% distinct() %>% rename('gene' = 'Linked_gene') %>% arrange(peak)
snp2 = snp2 %>% mutate(chr = sapply(peak, function(x) gsub('chr','', strsplit(peak,'-')[[1]][1]) %>% as.numeric))
snp2$start = sapply(snp2$peak, function(x){ strsplit(x,'-')[[1]][2] }) %>% as.numeric
snp2$end = sapply(snp2$peak, function(x){ strsplit(x,'-')[[1]][3] }) %>% as.numeric
snp2 = snp2 %>% filter(tf != 'NA') %>% arrange(chr, start,end)

## p2g

p2g = snp2 %>% dplyr::select(rsID,peak,Linked_gene) %>% distinct 
head(p2g,3)

p2g1 = snp2 %>% dplyr::select(rsID,peak) %>% distinct
p2g1$Linked_gene = NA
for (n in 1:nrow(p2g1)){
    gene = p2g[p2g$rsID == p2g1$rsID[n], ]$Linked_gene %>% paste0(., collapse = ',')
    p2g1[n, ]$Linked_gene = gene
}

## TF

svm = snp2 %>% dplyr::select(rsID,tf,preferred_allele) %>% distinct
svm = svm %>% pivot_wider(., values_from = tf, names_from = preferred_allele)
colnames(svm) = c('rsID','TF Loss','TF Gain')
head(svm,3)

result = left_join(p2g1, svm, by = 'rsID') %>% distinct
result

## snp anno

library(clusterProfiler)
library("org.Hs.eg.db")
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

head(snp,3)

peak = snp %>% filter(rsID %in% result$rsID)
peak = peak[, c('rsID','Chr','Start','Start')] %>% distinct
colnames(peak) = c('rsID','chr','start','end')
# peak$end = peak$start+1
head(peak,3)
dim(peak)

snp17.gr = GRanges(seqnames = peak$chr,
                   ranges = IRanges(peak$start, peak$end))

library(ChIPseeker)
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- getTagMatrix(snp17.gr, windows=promoter)

peakAnno <- annotatePeak(
    snp17.gr,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db")

peakAnno1 = peakAnno@anno@elementMetadata@listData %>% as.data.frame
peakAnno1

peak = peak %>% mutate(annotation = peakAnno1$annotation, SYMBOL = peakAnno1$SYMBOL)
peak

result = left_join(result, peak, by = 'rsID')
head(result,3)

result1 = result %>% mutate(Chr = gsub('chr','',chr) %>% as.numeric) %>% arrange(Chr,start) %>% dplyr::select(rsID,peak,annotation,SYMBOL,Linked_gene,'TF Loss','TF Gain')
result1

result1[,'TF Loss'] = lapply(result1[,'TF Loss'], function(x){
    ifelse(is.null(x), '-', paste(x, collapse = ','))
}) %>% unlist
result1[,'TF Gain'] = lapply(result1[,'TF Gain'], function(x){
    ifelse(is.null(x), '-', paste(x, collapse = ','))
}) %>% unlist
result1

result1$annotation[grep('Intron', result1$annotation)] = 'Intron'
result1

write.csv(result1, file = 'Figure5-sup_table.17snp.csv')

