#' @title mirProcessML
#'
#' @description  data pre-processing for miRNA-seq read count files specific for machine learning (ML). This function only runs under unix or unix-like operating systems. see \code{\link{mirProcess}}.
#' @param wd Working directory where all the read count \code{.txt} files are stored. Default is the current working directory.
#' @param setType Type of output data set, training or test set. Options are \code{"training"} and \code{"test"}. Default is \code{"training"}.
#' @details Make sure to follow the fie name naming convention for the read count files: \code{ID_database_targettype.txt}
#' @return Outputs a list with merged read counts from mutliple files, with annotation. No merging inlcuded: see \code{\link{mirProcess}}.
#' @import parallel
#' @examples
#' \dontrun{
#' readcountML <- mirProcessML()
#' }
#' @export
mirProcessML <- function(wd = getwd(), setType = "training"){
  # setting up parallel computing (using parallel package)
  n_cores <- detectCores() - 1
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl)) # close connection

  # locate to the working directory
  setwd(wd)

  # import the files
  rtpw <- readline("enter the root password: ")
  system("sudo -kS sh -c 'ls | grep .txt > filenames'", input = rtpw) # call system conmand to extract the txt file name into a temporary file
  inputDfm <- read.table(file = "filenames", stringsAsFactors = FALSE) # read the content of the
  system("sudo -kS rm filenames", input = rtpw) # call system command to remove the temporary fle
  colnames(inputDfm) <- "org.fileName"
  inputDfm$fileName <- sapply(inputDfm$org.fileName, function(x)unlist(strsplit(x, "\\."))[[1]], simplify = TRUE) # remove the extension of the file names
  inputDfm$targetType <- sapply(inputDfm$fileName, function(x)unlist(strsplit(x, "_"))[[3]], simplify =TRUE)
  inputDfm$targetType <- factor(inputDfm$targetType, levels = c(unique(inputDfm$targetType)))
  inputDfm$experimentID <- sapply(inputDfm$fileName, function(x)unlist(strsplit(x, "_"))[[1]], simplify = TRUE)

  if (setType == "test"){
    # parse the information and create a raw data list
    rawLstML <- parSapply(cl, inputDfm$fileName, function(x){
    temp <- read.table(file = paste(x, ".txt", sep=""), header = FALSE, stringsAsFactors = FALSE,
                         row.names = NULL)
    temp <- temp[-1,]
    colnames(temp)[1] <- "rawCount"
    colnames(temp)[2] <- unlist(strsplit(x, "_"))[2]

    row.names(temp) <- temp[,2]
    return(temp)
    }, simplify = FALSE, USE.NAMES = TRUE)
  }

  if (setType == "training"){
    rawLstML <- parSapply(cl, inputDfm$fileName, function(x){
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
    temp$miRNA_class <- factor(temp$miRNA_class, levels = c(unique(temp$miRNA_class)))
    return(temp)
    }, simplify = FALSE, USE.NAMES = TRUE)
  }

  return(rawLstML)
}

#' @title hairpinSet
#'
#' @description Produce a \code{.fasta} file as the training or test hairpin data set from the raw read count files. This function only works under Unix or Unix-like operating systems, with \code{pcregrep} installed.
#' @param dataLst Data list produced by \code{\link{mirProcessML}}.
#' @param setType The type of the output \code{fasta} file. Options are \code{"training"} or \code{"test"}. Default is \code{"training"}.
#' @return Outputs a \code{.fasta} file that can be used as the training set or test set for novel miRNA discovery. File names: training set \code{trainingSet.fasta}; test set \code{testSet.fasta}
#' @details Although the function is optimized for parallel computing, working with large reference file might stil result in long running time or system freezing, depending on the hardware configureation (CPU cores and RAM). Be sure to follow the file name naming convetion for refrence files: training reference \code{trainingRef.fasta}; test reference \code{testRef.fasta}.
#' @import parallel
#' @importFrom data.table rbindlist
#' @examples
#' \dontrun{
#' mylist <- mirProcessML() # produce data list
#' hairpinSet(mylist, setType = "training") # produce the hairpin training set
#' }
#' @export
hairpinSet <- function(dataLst, setType = "training"){
    # setting up parallel computing (using parallel package)
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl)) # close connection
    rtpw <- readline("enter the root password: ")

    if (setType == "training"){
      # prepare the hairpin vector from the raw data list
      tgtType <- unique(sapply(names(dataLst),
                               function(x)unlist(strsplit(x, "_"))[3], simplify = TRUE)) # extract the unique target types
      mergedLst <- parSapply(cl, tgtType, function(x){
      tempLst <- dataLst[grep(x, names(dataLst))]
      Dfm <- Reduce(function(i, j)merge(i[, c(2, 4)], j[, c(2, 4)], by = "miRNA_class", all = TRUE), tempLst)
      names(Dfm)[-1] <- sapply(strsplit(names(tempLst), "_"), "[[", 1) # use the function "[[" and the argument ", 1" to select the first element of the list element
      Dfm[is.na(Dfm) == TRUE] <- 0
      Dfm <- unique(Dfm)
      return(Dfm)
      }, simplify = FALSE, USE.NAMES = TRUE)

      dfm <- data.frame(mergedLst[[1]], stringsAsFactors = FALSE) # extract the hairpin lists and convert to dataframe
      dfm <- dfm[, -1]
      dfm <- as.matrix(dfm)
      V <- dfm[,1]
      V <- V[which(!V == 0)]
      V <- unique(V)
      write.table(V, file = "trainingIds.tmp", quote = FALSE, row.names = FALSE, col.names = FALSE)

      # output fastq file as the training set
      system("sudo -kS sh -c './training.smir'", input = rtpw)
      system("sudo -kS sh -c 'rm trainingIds.tmp'", input = rtpw)
    }

    if (setType == "test"){
      # prepare the hairpin vector from the raw data list
      dfm <- as.data.frame(rbindlist(dataLst), stringAsFactors = FALSE) # from package data.table
      V <- dfm[, 2]
      V <- V[which(!V == 0)]
      V <- unique(V)

      # output fastq file as the training set
      write.table(V, file = "testIds.tmp", quote = FALSE, row.names = FALSE, col.names = FALSE)

      # output fastq file as the training set
      system("sudo -kS sh -c './test.smir'", input = rtpw)
      system("sudo -kS sh -c 'rm testIds.tmp'", input = rtpw)
    }
}
