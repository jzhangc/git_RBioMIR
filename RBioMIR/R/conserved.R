.onAttach <- function(libname, pkgname) {
  packageStartupMessage("
                        Written by Jing Zhang, Ph.D. Please direct questions to jzhangcad@gmail.com.

                        (test only) For novel miRNA discovery, be sure to copy training.smir and test.smir files
                        to your working directory, and install pcregrep and parallel commands.")
  return(TRUE)
}

#' @title mirProcess
#'
#' @description  data pre-processing for miRNA-seq read count files. This function only runs under Unix or Unix-like operating systems. See \code{\link{mirProcessML}}.
#' @param wd Working directory where all the read count \code{.txt} files are stored. Default is the current working directory.
#' @details Make sure to follow the file name naming convention for the read count files: \code{ID_database_targettype.txt}
#' @return Outputs a list with merged read counts from mutliple files, with annotation.
#' @import parallel
#' @examples
#' \dontrun{
#' readcountMerged <- mirProcess()
#' }
#' @export
mirProcess <- function(wd = getwd()){
  # setting up parallel computing
  n_cores <- detectCores() - 1
  cl <- makeCluster(n_cores)

  # locate to the working directory
  setwd(wd)

  # import the files
  rtpw <- getPass("enter the root password: ")
  system("sudo -kS ls | grep .txt > filenames", input = rtpw) # call system conmand to extract the txt file name into a temporary file
  inputDfm <- read.table(file = "filenames", stringsAsFactors = FALSE) # read the content of the
  system("sudo -kS rm filenames", input = rtpw) # call system command to remove the temporary fle
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

  # combine the dataframes to generate count table for stats #
  # sidenote: seq(mylist) has the same effect of seq(length(mylist))
  # ave(): similar to tapply(), but WAY faster and doesn't have the length problem.
  maxdataLst <- sapply(names(rawdataLst),
                       function(x)rawdataLst[[x]][rawdataLst[[x]]$rawCount == ave(rawdataLst[[x]]$rawCount,
                                                                                  rawdataLst[[x]]$miRNA_class, FUN = max),],
                       simplify = FALSE, USE.NAMES = TRUE)

  tgtType <- unique(sapply(names(rawdataLst),
                           function(x)unlist(strsplit(x, "_"))[3], simplify = TRUE)) # extract the unique target types

  maxdataLstMg <- sapply(tgtType, function(x){
    tempLst <- lapply(maxdataLst[grep("miRNA", names(maxdataLst))], "[", c(1, 4)) # subset
    tempLst <- lapply(tempLst, unique)
    for (i in 1:length(tempLst)){
      names(tempLst[[i]])[1] <- unlist(strsplit(names(tempLst[i]), "_"))[1]
    }
    Dfm <- Reduce(function(i, j)merge(i, j, by = "miRNA_class", all = TRUE), tempLst)
    Dfm[is.na(Dfm) == TRUE] <- 0
    Dfm <- unique(Dfm)
    return(Dfm)
  }, simplify = FALSE, USE.NAMES = TRUE)

  maxdataLstMg <- maxdataLstMg[[1]]

  return(maxdataLstMg)
}


#' @title mirNrm
#'
#' @description  data Voom normalization for plotting
#' @param dfm Input dataframe.
#' @param wgt If or not to apply sample weight. Default is \code{FALSE}.
#' @return Outputs a dataframe with voom normalization
#' @importFrom limma voom voomWithQualityWeights
#' @importFrom edgeR DGEList calcNormFactors
#' @examples
#' \dontrun{
#' simNrm <- mirnrm(simDfm, wgt = TRUE)
#' }
#' @export
mirNrm <- function(dfm, wgt = FALSE){
  Nrm <- DGEList(counts = dfm[, -1], genes = dfm[1])
  Nrm <- calcNormFactors(Nrm)
  if (wgt == FALSE){
    Nrm <- voom(Nrm)
  } else {
    Nrm <- voomWithQualityWeights(Nrm)
  }
  NrmDfm <- data.frame(Nrm)
  return(NrmDfm)
}


#' @title mirFit
#'
#' @description Linear fitting and emperical Bayesian statistical test, with the capability of producing volcano distribution, among other useful plots.
#' @param dfm Input dataframe.
#' @param anno Annoation file for the samples. Format is \code{csv} and make sure the second column is for the conditions.
#' @param fileName output file name. Be sure to use quatation
#' @param wgt If or not to apply sample weight. Default is \code{FALSE}.
#' @return Outputs a list with limma eBayes fitting, QC plots and a volcano distribution
#' @importFrom limma voom voomWithQualityWeights plotMA plotSA lmFit eBayes topTable plotMDS
#' @importFrom edgeR DGEList calcNormFactors
#' @examples
#' \dontrun{
#' fitRNA <- mirFit(simDfm, wgt = TRUE)
#' }
#' @export
mirFit <- function(dfm, anno, fileName, wgt = FALSE){
  # load annoation
  SampleIndex <- read.csv(file = "condition.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  names(SampleIndex)[2] <- "Condition"

  # normalization
  Nrm <- DGEList(counts = dfm[, -1], genes = dfm[1])
  Nrm <- calcNormFactors(Nrm)
  plotMDS(Nrm)

  # voom normalization
  if (wgt == FALSE){
    Nrm <- voom(Nrm)
  } else {
    Nrm <- voomWithQualityWeights(Nrm)
  }

  # desgin matrix
  Exp<-factor(SampleIndex$Condition,levels=unique(SampleIndex$Condition))
  design<-model.matrix(~Exp)

  # eBayes fitting
  fitRNA<-lmFit(Nrm, design) # linear fitting
  fitRNA<-eBayes(fitRNA)
  plotSA(fitRNA,xlab="Average log-expression",ylab="log2(sigma)",
         zero.weights=FALSE, pch=16, cex=0.2) # Simga vs Average plot
  plotMA(fitRNA) # log-ratios vs mean average plot


  out <- list(fitted = fitRNA, Normalized_dataframe = data.frame(Nrm))

  exportTable <- topTable(fit = out[[1]], n = Inf, coef = attributes(design)[[2]][[2]][2], sort.by = "p", adjust.method = "fdr")

  # volcano plot
  with(exportTable, plot(logFC, -log10(adj.P.Val), pch = 20, cex = 0.5, main = "Volcano plot"))

  # Add colored points: red if padj < 0.05, orange of log2FC > 1, green if both)
  with(subset(exportTable, adj.P.Val < 0.05 ), points(logFC, -log10(adj.P.Val), pch = 20, cex = 0.5, col = "red"))
  with(subset(exportTable, abs(logFC) >= log2(1.5)), points(logFC, -log10(adj.P.Val), pch = 20, cex = 0.5,
                                                            col = "orange"))
  with(subset(exportTable, adj.P.Val<.05 & abs(logFC) >= log2(1.5)), points(logFC, -log10(adj.P.Val),
                                                                            pch = 20, cex = 0.5, col = "green"))

  # output
  write.csv(exportTable, file = paste(fileName, ".csv", sep = ""), row.names = FALSE)
  return(out)
}


#' @title mirPlot
#'
#' @description  Histogram. Also usable in most of the situations, expecially useful when end = 1.
#' @param dfm Input dataframe.
#' @param xTxtSize Font size of x-axis tick label. Default is 10.
#' @param yTxtSize Font size of y-axis tick label. Default is 10.
#' @param plotWidth Width of the output image file. Default is 170.
#' @param plotHeight Height of the output image file. Default is 150.
#' @return Histogram of the miRNA representative expression levels
#' @importFrom reshape2 melt
#' @importFrom grid grid.newpage grid.draw
#' @importFrom gtable gtable_add_cols gtable_add_grob
#' @import ggplot2
#' @examples
#' \dontrun{
#' mirPlot(miRNANrm, xTxtSize = 6, plotWidth = 360)
#' }
#' @export
mirPlot <- function(dfm, xTxtSize = 10, yTxtSize =10,
                    plotWidth = 170, plotHeight = 150){
  dfmplt <- melt(dfm, id.vars = colnames(dfm)[1]) # melt to make the dataframe for plotting

  # plot
  loclEnv <- environment()
  plt <- ggplot(dfmplt, aes(x = miRNA_class, y = value, fill= variable), environment = loclEnv) +
    geom_bar(position="dodge",stat="identity",color="black")+
    ggtitle(NULL) +
    xlab(NULL)+ # we can hide it using NULL
    ylab("Normalized read counts")+
    scale_y_continuous(expand = c(0, 0),
                       limits = c(0, with(dfmplt, max(value) * 1.1)))+
    theme(panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          legend.position = "bottom",
          legend.title = element_blank(),
          axis.text.x = element_text(size = xTxtSize, angle = 90, hjust = 1),
          axis.text.y = element_text(size = yTxtSize, hjust = 0.5)) +
    scale_fill_grey(start = 0)

  ## add the right-side y axis
  grid.newpage()

  # extract gtable
  pltgtb <- ggplot_gtable(ggplot_build(plt))

  # add the right side y axis
  Aa <- which(pltgtb$layout$name == "axis-l")
  pltgtb_a <- pltgtb$grobs[[Aa]]
  axs <- pltgtb_a$children[[2]]
  axs$widths <- rev(axs$widths)
  axs$grobs <- rev(axs$grobs)
  axs$grobs[[1]]$x <- axs$grobs[[1]]$x - unit(1, "npc") + unit(0.08, "cm")
  Ap <- c(subset(pltgtb$layout, name == "panel", select = t:r))
  pltgtb <- gtable_add_cols(pltgtb, pltgtb$widths[pltgtb$layout[Aa, ]$l], length(pltgtb$widths) - 1)
  pltgtb <- gtable_add_grob(pltgtb, axs, Ap$t, length(pltgtb$widths) - 1, Ap$b)

  # export the file and draw a preview
  ggsave(filename = paste(deparse(substitute(dfm)),".plot.pdf", sep = ""), plot = pltgtb,
         width = plotWidth, height = plotHeight, units = "mm",dpi = 600) # deparse(substitute(dfm)) converts object name into a character string
  grid.draw(pltgtb) # preview
}
