#' Save SCEmbroider files
#'
#'
#' Save a Seurat Object into SCEmbroider files. These file is used when loading new data in SCEmbroider program. Runs PCA, TSNE, UMAP if do not exists.
#' @param path Path of directory to save files
#' @param object A Seurat object
#' @param assay Name of Assay to make files
#' @param dims 	Which dimensions to use as input features (Only if tsne or umap is missing)
#' @param cluster.id If set, save name of the cluster idents with this column instead of default idents
#' @return A matrix of the infile
#' @export

saveSCEFile <- function(path ,object, assay = "RNA", dims = 1:5, cluster.id = NULL) {
  if(is.null(object@reductions[["pca"]])) {
    message("Cannot find 'pca' in this Seurat object\n Running RunPCA...")
    object <- RunPCA(object, assay = assay)
  }

  if(is.null(object@reductions[["tsne"]])) {
    message("Cannot find 'tsne' in this Seurat object\n Running RunTSNE...")
    object <- RunTSNE(object, assay = assay, dims = dims)
  }

  if(is.null(object@reductions[["umap"]])) {
    message("Cannot find 'umap' in this Seurat object\n Running RunUMAP...")
    object <- RunUMAP(object, assay = assay, dims = dims)
  }

  TSNE <- object@reductions[["tsne"]]@cell.embeddings
  UMAP <- object@reductions[["umap"]]@cell.embeddings


  GCidx <- object@assays[[assay]]@data@i
  GCnum <- object@assays[[assay]]@data@p
  GCval <- object@assays[[assay]]@data@x
  GCdim <- object@assays[[assay]]@data@Dim

  nFeature_name = paste("nFeature", assay, sep = "_")
  GCsum <- object@meta.data[[nFeature_name]]
  GENE  <- object@assays[[assay]]@data@Dimnames[[1]]
  BCODE <- object@assays[[assay]]@data@Dimnames[[2]]
  META <- object[[]]
  Annotation <- META %>% dplyr::select(orig.ident)

  if(!is.null(cluster.id)) {
    Annotation$cluster.ident <- META$cluster.id
  }
  else {
    Annotation$cluster.ident <- Idents(object = object)
  }
  writeMat('./SCE/GC_Data.mat',
           GCidx = GCidx, GCval = GCval, GCnum = GCnum, GCdim = GCdim, GCsum = GCsum, TSNE = TSNE, UMAP = UMAP)

  write.csv(GENE, './SCE/Gene_Name.csv')
  write.csv(Annotation, './SCE/Cell_Type.csv')
}
