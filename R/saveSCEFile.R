#' Save SCEmbroider files
#'
#'
#' Save a Seurat Object into SCEmbroider files. These file is used when loading new data in SCEmbroider program. Runs PCA, TSNE, UMAP if do not exists.
#'
#' @param object A Seurat object
#' @param path Path of directory to save files
#' @param assay Name of Assay to make files
#' @param dims 	Which dimensions to use as input features (Only if pca, tsne or umap is missing)
#' @param cluster.id If set, save name of the cluster idents with this column instead of default idents
#' @return A Seurat object modified in this function
#' @export

saveSCEFile <- function(object, path, assay = "RNA", dims = 1:5, cluster.id = NULL) {
  new_object <- UpdateSeuratObject(object)

  if(length(VariableFeatures(new_object)) == 0) {
    message("\n\nCannot find VariableFeatures in this Seurat object\n Running SCTransform...")
    new_object <- SCTransform(new_object, assay = assay, verbose = TRUE)
    assay = "SCT"
  }

  if(is.null(new_object@reductions[["pca"]])) {
    message("\n\nCannot find 'pca' in this Seurat object\n Running RunPCA...")
    new_object <- RunPCA(new_object, verbose = FALSE)
  }

  META <- new_object[[]]
  if(is.null(cluster.id) & is.null(META$seurat_clusters)){
    yesNo <- askYesNo("\n\nCannot find cluster info in this object. Do you want to find
                      clusters automatically?")
    if(yesNo == TRUE) {
      new_object <- FindNeighbors(new_object, reduction = "pca", dims = dims)
      new_object <- FindClusters(new_object, verbose = FALSE)
    } else {
      return(new_object)
    }
  }

  if(is.null(new_object@reductions[["tsne"]])) {
    message("\n\nCannot find 'tsne' in this Seurat object\n Running RunTSNE...")
    new_object <- RunTSNE(new_object, assay = assay, dims = dims)
  }

  if(is.null(new_object@reductions[["umap"]])) {
    message("\n\nCannot find 'umap' in this Seurat object\n Running RunUMAP...")
    new_object <- RunUMAP(new_object, assay = assay, dims = dims)
  }

  TSNE <- new_object@reductions[["tsne"]]@cell.embeddings
  UMAP <- new_object@reductions[["umap"]]@cell.embeddings


  GCidx <- new_object@assays[[assay]]@data@i
  GCnum <- new_object@assays[[assay]]@data@p
  GCval <- new_object@assays[[assay]]@data@x
  GCdim <- new_object@assays[[assay]]@data@Dim

  nFeature_name = paste("nFeature", assay, sep = "_")
  GCsum <- new_object@meta.data[[nFeature_name]]
  GENE  <- new_object@assays[[assay]]@data@Dimnames[[1]]
  BCODE <- new_object@assays[[assay]]@data@Dimnames[[2]]
  META <- new_object[[]]
  Annotation <- META %>% dplyr::select(orig.ident)

  if(!is.null(cluster.id)) {
    Annotation$cluster.ident <- dplyr::select(META, cluster.id)
  }
  else {
    Annotation$cluster.ident <- Idents(object = new_object)
  }
  writeMat('./SCE/GC_Data.mat',
           GCidx = GCidx, GCval = GCval, GCnum = GCnum, GCdim = GCdim, GCsum = GCsum, TSNE = TSNE, UMAP = UMAP)

  write.csv(GENE, './SCE/Gene_Name.csv')
  write.csv(Annotation, './SCE/Cell_Type.csv')

  return(new_object)
}

