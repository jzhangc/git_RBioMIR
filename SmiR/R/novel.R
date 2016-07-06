#' @title mirProcessML
#'
#' @description  data pre-processing for miRNA-seq read count files specific for machine learning (ML). This function only runs under unix or unix-like operating systems. see \code{\link{mirProcess}}.
#' @param wd Working directory where all the read count \code{.txt} files are stored. Default is the current working directory.
#' @details Make sure to follow the fie name naming convention for the read count files: ID_database_targettype.txt
#' @return Outputs a list with merged read counts from mutliple files, with annotation. No merging inlcuded: see \code{\link{mirProcess}}.
#' @examples
#' \dontrun{
#' readcountML <- mirProcessML()
#' }
#' @export
mirProcessML <- function(wd = getwd()){
  # locate to the working directory
  setwd(wd)

  # import the files
  system("sudo -kS sh -c 'ls | grep .txt > filenames'", input = readline("Enter your password: ")) # call system conmand to extract the txt file name into a temporary file
  inputDfm <- read.table(file = "filenames", stringsAsFactors = FALSE) # read the content of the
  system("sudo -kS rm filenames", input = readline("Enter your password: ")) # call system command to remove the temporary fle
  colnames(inputDfm) <- "org.fileName"
  inputDfm$fileName <- sapply(inputDfm$org.fileName, function(x)unlist(strsplit(x, "\\."))[[1]], simplify = TRUE) # remove the extension of the file names
  inputDfm$targetType <- sapply(inputDfm$fileName, function(x)unlist(strsplit(x, "_"))[[3]], simplify =TRUE)
  inputDfm$targetType <- factor(inputDfm$targetType, levels = c(unique(inputDfm$targetType)))
  inputDfm$experimentID <- sapply(inputDfm$fileName, function(x)unlist(strsplit(x, "_"))[[1]], simplify = TRUE)

  # parse the information and create a raw data list
  rawdataLst <- sapply(inputDfm$fileName, function(x){
    temp <- read.table(file = paste(x, ".txt", sep=""), header = FALSE, stringsAsFactors = FALSE,
                       row.names = NULL)
    temp <- temp[-1,]
    colnames(temp)[1] <- "rawCount"
    colnames(temp)[2] <- unlist(strsplit(x, "_"))[3]

    row.names(temp) <- temp[,2]

    temp$species <- sapply(temp[[2]], function(i)unlist(strsplit(i, "-"))[1], simplify = TRUE)

    temp$miRNA_class <- sapply(temp[[2]], function(i)paste(unlist(strsplit(i, "-"))[2],
                                                           "-",
                                                           unlist(strsplit(i, "-"))[3],
                                                           sep = ""),
                               simplify = TRUE)

    temp$species <- factor(temp$species, levels = c(unique(temp$species)))
    temp$miRNA_class <- gsub("-NA", "", temp$miRNA_class) # remove the -NA string that appears when the mirbase id is like xxx-miRXXX
    temp$miRNA_class <- factor(temp$miRNA_class, levels = c(unique(temp$miRNA_class)))
    return(temp)
  }, simplify = FALSE, USE.NAMES = TRUE)

  return(rawdataLst)
}

#' @title hairpinSets
#'
#' @description Produce a \code{.fasta} file as the training hairpin data set from the raw
#' @param dataLst If or not to apply sample weight.
#' @param refFileName Full directory and file name of ther reference.
#' @return Outputs a \code{.fasta} file that can be used as the training set for novel miRNA discovery.
#' @details Working with large reference file might result in long running time or system freezing, depending on the hardware configureation (mainly RAM). It is recommanded to build an index for the reference file prior to this operation when the file is large (multi-GB).
#' @importFrom seqinr read.fasta
#' @examples
#' \dontrun{
#' hairpinTrainset(refFileName = "~/OneDrive/Storey lab/current_work/miRNA_pipeline/reference/hairpin.fa", rawdataLst) # produce the hairpin training set
#' }
#' @export
hairpinSets <- function(refFileName, dataLst){
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
