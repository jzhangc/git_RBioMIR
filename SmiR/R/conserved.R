.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Written by Jing Zhang, Ph.D. Please direct questions to jzhangcad@gmail.com.")
  return(TRUE)
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
#' simNrm <- mirnrm(simDfm, wgt = TRUE)#'
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
#' @description  data Voom normalization for plotting
#' @param dfm Input dataframe.
#' @param fileName output file name. Be sure to use quatation
#' @param wgt If or not to apply sample weight. Default is \code{FALSE}.
#' @return Outputs a list with limma eBayes fitting, QC plots and a volcano distribution
#' @importFrom limma voom voomWithQualityWeights plotMA plotSA lmFit eBayes topTable
#' @importFrom edgeR DGEList calcNormFactors
#' @examples
#' \dontrun{
#' fitRNA <- mirFit(simDfm, wgt = TRUE)#'
#' }
#' @export
mirFit <- function(dfm, fileName, wgt = FALSE){
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
  SampleIndex <- data.frame(Condition = sapply(colnames(Nrm), function(x)substr(x, nchar(x)-1, nchar(x))),
                            Sample = colnames(Nrm), stringsAsFactors = FALSE)

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
