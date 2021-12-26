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
#' @slot sample_groups factor. Primary sample group factor.
#' @slot miRNA_database str. MiR database.
#' @slot selected_species str. The miRBase species code for the selected species of interest.
#' @slot total_species str. Total number of species (in miRBase species code).
#' @slot files_processed str. Processed file(s) name containing raw reads.
#' @slot target_annotation_file_processed str. Processed sample annotation file name.
#'
#' @details 1. The class links the column of \code{raw_read_count}, rows of \code{targets}, \code{sample_library_sizes}, \code{sample_groups} and \code{genes}.
#'
#'          2. Subsetting in "1" does not subset \code{genes_complete_annotation} as this data.frame may contain replicates.
mircount <- setClass("mircount", representation(
  raw_read_count = "dfORmatrix",
  sample_library_sizes = "numeric",
  genes_complete_annotation = "data.frame",
  working_gene_annot_var_name = "character",
  genes = "character",
  targets = "data.frame",
  sample_groups_var_name = "character",
  sample_groups = "factor",
  miRNA_database = "character",
  selected_species = "character",
  total_species = "character",
  files_processed = "character",
  target_annotation_file_processed = "character"
))

#' @title dim.mircount
#' @description dim method for \code{mircount} class.
#'
#' @details dim is the dim of \code{raw_read_count}.
#' @export
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
#' @export
length.mircount <- function(x) {
  dim(x@raw_read_count)[2]
}

#' @title as.matrix.mircount
#' @description as.matrix method for \code{mircount} class.
#'
#' @details It converts \code{raw_read_count} into a matrix
#' @export
as.matrix.mircount <- function(x) {
  y <- as.matrix(x@raw_read_count)
  return(y)
}

#' @title as.data.frame.mircount
#' @description as.data.frame method for \code{mircount} class.
#'
#' @details It converts \code{raw_read_count} into a data.frame.
#' @export
as.data.frame.mircount <- function(x) {
  y <- as.data.frame(x@raw_read_count)
  return(y)
}

#' @title \code{mircount} substting method
#' @description Automatic subsetting method for \code{mircount} class.
#'
#' @details It links the column of \code{raw_read_count}, rows of \code{targets}, \code{sample_library_sizes} and \code{genes}.
#' @examples
#' \dontrun{
#' a <- new("mircount")
#' b <- a[1:5, 1:5]  # leads to raw_read_count[1:5, 1:5] and targets[1:5, ]
#' }
#' @export
setMethod("[", "mircount", definition = function(x, i, j, ..., drop = TRUE){
  initialize(x,
             raw_read_count = x@raw_read_count[i, j, ..., drop = drop],
             genes = x@genes[i],
             sample_library_sizes = x@sample_library_sizes[j],
             sample_groups = factor(x@sample_groups[j], levels = unique(x@sample_groups[j])),
             targets = x@targets[j, , ..., drop = drop])
})

#' @title as.mircount
#' @description Convert compatible S3 classes to an \code{mircount} class.
#' @export
as.mircount <- function(object, ...) {
  UseMethod("as.mircount", object)
}

#' @title as.mircount.mir_count
#' @rdname as.mircount
#' @description \code{as.mircount} for S3 \code{mir_count} method
#' @method as.mircount mir_count
#' @export
as.mircount.mir_count <- function(x) {
  y <- as.mircount.default(raw_read_count = x$raw_read_count,
                           sample_library_sizes = x$sample_library_sizes,
                           genes_complete_annotation = x$genes_complete_annotation,
                           working_gene_annot_var_name = x$working_gene_annot_var_name,
                           genes = x$genes,
                           targets = x$targets,
                           sample_groups_var_name = x$sample_groups_var_name,
                           sample_groups = x$sample_groups,
                           miRNA_database = x$miRNA_database,
                           selected_species = x$selected_species,
                           total_species = x$total_species,
                           files_processed = x$files_processed,
                           target_annotation_file_processed = x$target_annotation_file_processed)
  return(y)
}

#' @title as.mircount.default
#' @rdname as.mircount
#' @description Default method for \code{as.mircount}
#' @method as.mircount default
#' @export
as.mircount.default <- function(raw_read_count, sample_library_sizes,
                                genes_complete_annotation, working_gene_annot_var_name,
                                genes, targets, sample_groups_var_name,
                                sample_groups, miRNA_database, selected_species,
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
  y@sample_groups = sample_groups
  y@miRNA_database = miRNA_database
  y@selected_species = selected_species
  y@total_species = total_species
  y@files_processed = files_processed
  y@target_annotation_file_processed = target_annotation_file_processed
  return(y) # return the mapped y
}

#' @export
print.mircount <- function(x, ...){
  cat("MiRNA raw reads processing summary:\n")
  cat("\n")
  cat(paste0(" MiRNA database: ", x@miRNA_database, "\n"))
  cat(paste0(" Selected species: ", paste0(x@selected_species, collapse = " "), "\n"))
  cat(paste0(" Total number ", x@working_gene_annot_var_name, " targets (no merging): ", nrow(x@genes_complete_annotation), "\n"))
  cat(paste0(" Total number of unique ", x@working_gene_annot_var_name, " targets ", length(x@genes), "\n"))
  cat("\n")
  cat(paste0(" Files read: ", "\n"))
  cat(paste0(" ", x@files_processed, collapse = "\n"))
  cat("\n\n")
}
