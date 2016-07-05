#' @title hairpinTraining
#'
#' @description Produce a \code{.fasta} file as the training hairpin data set from the raw
#' @param refFileName File name of ther reference.
#' @param wgt If or not to apply sample weight. Default is \code{FALSE}.
#' @return Outputs a \code{.fasta} file that can be used as the training set for novel miRNA discovery.
#' @importFrom seqinr read.fasta
#' @examples
#' \dontrun{
#' hairpinTrainset(refFileName = "~/OneDrive/Storey lab/current_work/miRNA_pipeline/reference/hairpin.fa", rawdataLst) # produce the hairpin training set
#' }
#' @export
hairpinTraining <- function(refFileName, dataLst){
  # load the original miRBase hairpin fasta files (reqinr::read.fasta) and convert RNA sequences into DNA sequences (U > T)
  ref <- read.fasta(file = refFileName, as.string = TRUE, forceDNAtolower = FALSE)

  ref <- lapply(seq(ref), function(x)gsub("U", "T", ref[[x]])) # replace "U" with "T"

  # extract attributes from the list
  AttN <- sapply(seq(ref), function(i)attributes(ref[[i]])[[1]])
  AttAno <- sapply(seq(ref), function(i)attributes(ref[[i]])[[2]])
  refSeq <- unlist(ref)

  refDfm <- data.frame(hairpin = AttN, sequence = refSeq, annot = AttAno, stringsAsFactors = FALSE) # replace "hairpin"

  # prepare the hairpin vector from the raw data list
  tgtType <- unique(sapply(names(dataLst),
                           function(x)unlist(strsplit(x, "_"))[3], simplify = TRUE)) # extract the unique target types
  mergedLst <- sapply(tgtType, function(x){
    tempLst <- dataLst[grep(x, names(dataLst))]
    Dfm <- Reduce(function(i, j)merge(i[, c(2, 4)], j[, c(2, 4)], by = "miRNA_class", all = TRUE), tempLst)
    names(Dfm)[-1] <- sapply(strsplit(names(tempLst), "_"), "[[", 1) # use the function "[[" and the argument ", 1" to select the first element of the list element
    Dfm[is.na(Dfm) == TRUE] <- 0
    Dfm <- unique(Dfm)
    return(Dfm)
  }, simplify = FALSE, USE.NAMES = TRUE)

  dfm <- data.frame(mergedLst[[1]], stringsAsFactors = FALSE) # extract the hairpin lists and convert to dataframe
  dfm <- dfm [, -1]
  dfm <- as.matrix(dfm)
  V <- dfm[,1]
  V <- V[which(!V == 0)]
  V <- unique(V)

  # output fastq file as the training set
  out <- subset(refDfm, hairpin %in% V)

  out <- sapply(seq(nrow(out)), function(x){
    tmpV <- rbind(out$annot[[x]], out$sequence[[x]])
    return(tmpV)
  }
  )
  out <- as.vector(out)

  write.table(out, file = paste(deparse(substitute(dataLst)),".fastq", sep = ""),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
}
