#' class union for S4
dfORmatrix <- setClassUnion("dfORmatrix", c("data.frame", "matrix"))  # use more than one class

#' @title mircount
#' @description  S4 class equivalent to the S3 class \code{mir_count}.
#'
#' @slot raw_read_count data.frame or matrix. Raw miRNA read counts, with row: genes, column: samples.
#' @slot sample_library_sizes numeric. Vector for sample library sizes.
#' @slot genes_complete_annotation data.frame. Complete gene annotation.
#' @slot genes str. Vector for gene names.
#' @slot targets data.frame. Complete sample annotation data frame.
#' @slot sample_groups_var_name str. Variable name for sample grouping in the sample annotation.
#' @slot total_species str. Total number of species (in miRBase species code).
#' @slot files_processed str. Processed file(s) name containing raw reads.
#' @slot target_annotation_file_processed str. Processed sample annotation file name.
mircount <- setClass("mircount", representation(
  raw_read_count="dfORmatrix",
  sample_library_sizes="numeric",
  genes_complete_annotation="data.frame",
  working_gene_annot_var_name="character",
  genes="character",
  targets="data.frame",
  sample_groups_var_name="character",
  total_species="character",
  files_processed="character",
  target_annotation_file_processed="character"
))

#' @title dim.mircount
#' @description dim method for \code{mircount} class.
#'
#' @details dim is the dim of \code{raw_read_count}.
dim.mircount <- function(x) {
  if (is.null(x@raw_read_count)) {
    c(0, 0)
  } else {
    dim(x@raw_read_count)
  }
}

#' @title length.mircount
#' @description length method for \code{mircount} class.
#'
#' @details length is the number of samples
length.mircount <- function(x) {
  dim(x@raw_read_count)[2]
}

#' @title as.matrix.mircount
#' @description as.matrix method for \code{mircount} class.
#'
#' @details It converts \code{raw_read_count} into a matrix
as.matrix.mircount <- function(x) {
  y <- as.matrix(x@raw_read_count)
  return(y)
}

#' @title as.data.frame.mircount
#' @description as.data.frame method for \code{mircount} class.
#'
#' @details It converts \code{raw_read_count} into a data.frame.
as.data.frame.mircount <- function(x) {
  y <- as.data.frame(x@raw_read_count)
  return(y)
}

#' @title \code{mircount} substting method
#' @description Automatic subsetting method for \code{mircount} class.
#'
#' @details It links the column of \code{raw_read_count} and rows of \code{targets}.
#' @examples
#' \dontrun{
#' a <- new("mircount")
#' b <- a[1:5, 1:5]  # leads to raw_read_count[1:5, 1:5] and targets[1:5, ]
#' }
setMethod("[", "mircount", definition = function(x, i, j, ..., drop=TRUE){
  initialize(x, raw_read_count = x@raw_read_count[i, j, ..., drop=drop], targets = x@targets[j, , ..., drop=drop])
})

#' @title as.mircount
#' @description Convert compatible S3 classes to an \code{mircount} class.
#' @export
as.mircount <- function(object, ...) {
  UseMethod("as.mircount", object)
}

#' @title as.mircount.default
#' @description Default method for \code{as.mircount}
as.mircount.default <- function(raw_read_count, sample_library_sizes,
                                   genes_complete_annotation, working_gene_annot_var_name,
                                   genes, targets, sample_groups_var_name,
                                   total_species, files_processed, target_annotation_file_processed) {
  # covert the existing mir_count objects into the new S4 class
  y <- new("mircount")  # set up a new instance (y) fo the target class (i.e. "mircount")
  y@raw_read_count = raw_read_count
  y@sample_library_sizes = sample_library_sizes
  y@genes_complete_annotation = genes_complete_annotation
  y@working_gene_annot_var_name = working_gene_annot_var_name
  y@genes = genes
  y@targets = targets
  y@sample_groups_var_name = sample_groups_var_name
  y@total_species = total_species
  y@files_processed = files_processed
  y@target_annotation_file_processed = target_annotation_file_processed
  return(y) # return the mapped y
}

#' @title as.mircount.mir_count
#' @description \code{as.mircount} for S3 \code{mir_count} method for
as.mircount.mir_count <- function(x) {
  y <- as.mircount.default(raw_read_count=x$raw_read_count,
                              sample_library_sizes=x$sample_library_sizes,
                              genes_complete_annotation=x$genes_complete_annotation,
                              working_gene_annot_var_name=x$working_gene_annot_var_name,
                              genes=x$genes,
                              targets=x$targets,
                              sample_groups_var_name=x$sample_groups_var_name,
                              total_species=x$total_species,
                              files_processed=x$files_processed,
                              target_annotation_file_processed=x$target_annotation_file_processed)
  return(y)
}
