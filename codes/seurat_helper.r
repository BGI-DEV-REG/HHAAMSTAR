# Utilize the function in seurat for data analysis
library(Seurat)
library(SoupX)
library(scDblFinder)
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

# remove Ambient RNAs by contamination
runSoupX_Droplets <- function(matrix_path=NULL, expression_matrix=NULL, droplet_matrix=NULL){
    # loading data
    if (!is.null(matrix_path)){
        toc = Seurat::Read10X(file.path(matrix_path, "04.Matrix", "FilterMatrix"), gene.column = 1)
        tod = Seurat::Read10X(file.path(matrix_path, "02.cDNAAnno", "RawMatrix"), gene.column = 1)
    }
    else{
        toc=expression_matrix
        tod=droplet_matrix
    }
    #
    genes=intersect(rownames(toc),rownames(tod))
    toc=toc[genes,]
    tod=tod[genes,]
    sc = SoupChannel(tod, toc)
    # get clusters info
    obj <- CreateSeuratObject(counts = toc)
    obj=Process_RNA(obj)
    soupx_groups = Idents(obj)
    #
    sc = setClusters(sc, soupx_groups)
    sc = autoEstCont(sc, doPlot=FALSE)
    out = adjustCounts(sc,roundToInt = TRUE)
    return (out)
}

# export 10X format matrix of single-cell data
ExportData_10X_format <- function(data=NULL, out_dir=NULL, assay='RNA', if_compress=TRUE){
  #
  if (is.null(data)){stop('please supply the data')}
  #
  if(assay=='ATAC'){
    if (!dir.exists(out_dir)) {dir.create(out_dir)}
    writeLines(colnames(data), paste0(out_dir,'/barcodes.tsv'))
    writeLines(rownames(data), paste0(out_dir,'/peaks.bed'))
    Matrix::writeMM(data, file=paste0(out_dir,'/matrix.mtx'))
    if (if_compress){
        R.utils::gzip(paste0(out_dir,'/barcodes.tsv'))
        R.utils::gzip(paste0(out_dir,'/peaks.bed'))
        R.utils::gzip(paste0(out_dir,'/matrix.mtx'))
    }
  }
  #
  if(assay=='RNA'){
    if (!dir.exists(out_dir)) {dir.create(out_dir)}
    writeLines(colnames(data), paste0(out_dir,'/barcodes.tsv'))
    writeLines(rownames(data), paste0(out_dir,'/features.tsv'))
    Matrix::writeMM(data, file=paste0(out_dir,'/matrix.mtx'))
    if (if_compress){
        R.utils::gzip(paste0(out_dir,'/barcodes.tsv'))
        R.utils::gzip(paste0(out_dir,'/features.tsv'))
        R.utils::gzip(paste0(out_dir,'/matrix.mtx'))
    }
  }
  #
  print('Success to output data')
}

# get overlap genes from multiple seurat objects
Overlap_Seurat_Genes <- function(list_seurat=NULL){
    if (!inherits(x = list_seurat, what = "list")) {
        cli_abort(message = "{.code list_seurat} must be environmental variable of class {.val list}")
    }
    for (i in 1:length(x = list_seurat)) {
        if (!inherits(x = list_seurat[[i]], what = "Seurat")) {
            cli_abort("One or more of entries in {.code list_seurat} are not objects of class {.val Seurat}")
        }
    }
    gene_list <- lapply(X = list_seurat, FUN = function(x) {
        x <- rownames(x)
    })
    overlap_genes <- purrr::reduce(gene_list, function(x, y) {
        intersect(x, y)
    })
}

 Compute QC index of seurat object of RNA
QC_RNA <- function(proj=proj, genome='hg38', assay='RNA'){
  DefaultAssay(proj) <- assay
  #
  if (genome=='hg38'){
    proj[["percent.mt"]] <- PercentageFeatureSet(proj, pattern = "^MT-")
    proj[["percent.ribo"]] <- PercentageFeatureSet(proj, pattern = "^RP[SL]")
  }
  if (genome=='mm10'){
    features=grep(pattern = "^MT-", x = rownames(proj), value = TRUE,ignore.case = TRUE)
    proj[["percent.mt"]] <- PercentageFeatureSet(proj, pattern = "^mt-",features=features)
    proj[["percent.ribo"]] <- PercentageFeatureSet(proj, pattern = "^Rp[sl]")
  }
  #
  return(proj)
}

loadRDSData <- function(dir, sample_names, merge=TRUE, overlap_genes=TRUE){
    print (sample_names)
    file_paths=paste0(dir, sample_names, '.rds')
    for (i in file_paths){print (i); print (file.exists(i))}
    #
    seurat.list <- lapply(X = file_paths, FUN = function(x) {x=readRDS(x)})
    #
    if(overlap_genes){
        # get overlap genes
        genes=Overlap_Seurat_Genes(seurat.list)
        # extract genes of each seurat object
        seurat.list <- lapply(X = seurat.list, FUN = function(x) {x <- x[genes, ]})
    }
    #
    if(merge){
        #seurat.obj <- Merge_Seurat_List(list_seurat = seurat.list)
        seurat.obj <- merge(x = seurat.list[[1]], y = seurat.list[c(1:length(seurat.list))[-1]])
        #seurat.obj <- merge(x = seurat.list[[1]], y = seurat.list[c(1:length(seurat.list))[-1]], add.cell.ids = sample_names)
        return(seurat.obj)
    }
    else{
        return(seurat.list)
    }
}


