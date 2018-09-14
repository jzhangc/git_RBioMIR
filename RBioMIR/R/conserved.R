#' @title mirProcess
#'
#' @description Data pre-processing for miRNA-seq read count files. Also eee \code{\link{mirProcessML}}.
#' @param path Path to raw files. Default is the system working directory.
#' @param raw.file.sep Raw read count file separators. Default is \code{""\"\"}, i.e. white space.
#' @param species Species code, following the three letter naming convention.
#' @param target.annot.file Annotation file describing filenames and targets, and should be in \code{csv} format.
#' @param database MiRNA database, only for miRNA naming conventions. Currently the function only takes "mirbase".
#' @param parallelComputing Wether to use parallel computing or not. Default is \code{TRUE}.
#' @param cluterType clusterType Only set when \code{parallelComputing = TRUE}, the type for parallel cluster. Options are \code{"PSOCK"} (all operating systems) and \code{"FORK"} (macOS and Unix-like system only). Default is \code{"PSOCK"}.
#' @details When \code{species} left as \code{NULL}, the function will use all the species detected in the raw read count files.
#'          The raw count files are usually in \code{.txt} format with "read count" and "mirna" columns. The read count files can be obtained using the shell program \code{mirna_processing}.
#' @return Outputs a \code{mirna_count} with merged read counts from mutliple files, with annotation.
#' @import foreach
#' @import doParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @examples
#' \dontrun{
#' readcountMerged <- mirProcess()
#' }
#' @export
mirProcess <- function(path = getwd(), species = NULL, target.annot.file = NULL, database = "mirbase", raw.file.sep = "",
                       parallelComputing = FALSE, clusterType = "FORK"){
  ## check argument
  annot_name_length <- length(unlist(strsplit(target.annot.file, "\\.")))
  annot_ext <- unlist(strsplit(target.annot.file, "\\."))[annot_name_length]
  if (!database %in% c("mirbase")) stop("For now, the function only accepts database = \"mirbase\".")

  ## load files
  # check and load annotation
  if (is.null(target.annot.file) | annot_ext != "csv"){
    tgt <- NULL
  } else {
    cat("Loading target annotation file...")
    tgt <- read.csv(file = target.annot.file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    cat("Done!\n")
  }

  # set read files
  filename <- list.files(path = path, pattern = ".txt")
  filename_wo_ext <- sub("[.][^.]*$", "", filename)  # general expression to remove extension, i.e. a.b.c becomes a.b

  # load reads
  cat("Processing read count files...")
  raw_list <- vector(mode = "list", length = length(filename))
  if (!parallelComputing){ # single core
    raw_list[] <- foreach(i = filename) %do% {
      tmp <- read.table(file = paste0(path, "/", i), header = FALSE, sep = raw.file.sep, stringsAsFactors = FALSE,
                        col.names = c("read_count", "mirna"), row.names = NULL)
      tmp <- tmp[, c(2, 1)]

      tmp$species <- sapply(tmp$mirna, function(m)unlist(strsplit(m, "-"))[1], simplify = TRUE)
      tmp$mirna_class <- sapply(tmp$mirna, function(n)paste(unlist(strsplit(n, "-"))[2],
                                                            "-",
                                                            unlist(strsplit(n, "-"))[3],
                                                            sep = ""), simplify = TRUE)
      tmp
    }
  } else {  # parallel
    # set clusters
    n_cores <- detectCores() - 1
    cl <- makeCluster(n_cores, clusterType = clusterType)
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) # close connect when exiting the function

    # file processing
    raw_list[] <- foreach(i = filename, .packages = "foreach") %dopar% {
      tmp <- read.table(file = paste0(path, "/", i), header = FALSE, sep = raw.file.sep, stringsAsFactors = FALSE,
                        col.names = c("read_count", "mirna"), row.names = NULL)
      tmp <- tmp[, c(2, 1)]

      tmp$species <- sapply(tmp$mirna, function(m)unlist(strsplit(m, "-"))[1], simplify = TRUE)
      tmp$mirna_class <- sapply(tmp$mirna, function(n)paste(unlist(strsplit(n, "-"))[2],
                                                            "-",
                                                            unlist(strsplit(n, "-"))[3],
                                                            sep = ""), simplify = TRUE)
      tmp
    }
  }
  names(raw_list) <- filename_wo_ext
  cat("Done!\n")

  # file loaded message
  cat("\n")
  cat("Files loaded: \n")
  for (i in filename){
    cat(paste0("\t", i, "\n"))
  }

  # get and check species
  tot_species <- foreach(i = raw_list, .combine = "c") %do% i$species
  tot_species <- unique(tot_species)
  if (is.null(species) | !all(species %in% tot_species)){
    warning("Set species not all found in data, or species not set. Proceed with all species (may be very slow).")
    species <- tot_species
  }

  ## merge reads
  species_list <- vector(mode = "list", length = length(filename))
  species_list <- foreach(i = 1:length(filename)) %do% {
    tmpdfm <- raw_list[[i]][raw_list[[i]]$species %in% species, ]
    tmpdfm <- tmpdfm[, c("mirna", "read_count")]
    names(tmpdfm)[2] <- filename_wo_ext[i]
    tmpdfm
  }
  names(species_list) <- filename_wo_ext

  ## output
  out_dfm <- Reduce(function(i, j)merge(i, j, all = TRUE), species_list)
  out_dfm[is.na(out_dfm) == TRUE] <- 0
  out <- list(raw_read_count = out_dfm,
              target = tgt,
              selected_species = species,
              tot_species = tot_species)
  class(out) <- "mir_count"
  return(out)
}


#' @title mirNrm
#'
#' @description  data Voom normalization for plotting
#' @param dfm Input dataframe.
#' @param count_threshold Read count threshold. No filtering will be applied when set \code{"none"}. Otherwise, a numeric number can be set as the minimum read count for filtering. DDefault is \code{"none"}.
#' @param wgt If or not to apply sample weight. Default is \code{FALSE}.
#' @return Outputs a dataframe with voom normalization
#' @importFrom limma voom voomWithQualityWeights
#' @import edgeR
#' @examples
#' \dontrun{
#' simNrm <- mirnrm(simDfm, wgt = TRUE)
#' }
#' @export
mirNrm <- function(dfm, count_threshold = "none", wgt = FALSE){
  Nrm <- DGEList(counts = dfm[, -1], genes = dfm[1])
  # filtering
  if (count_threshold != "none"){ # set the count threshold for filtering
    count_s <- rowSums(Nrm$counts) # thresholdd
    isexpr <- count_s > count_threshold

    Nrm <- Nrm[isexpr, , keep.lib.size = FALSE] # filtering
  }

  # normalization
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
#' @param fileName output file name. Be sure to use quatation.
#' @param count_threshold Read count threshold. No filtering will be applied when set \code{"none"}. Otherwise, a numeric number can be set as the minimum read count for filtering. DDefault is \code{"none"}.
#' @param wgt If or not to apply sample weight. Default is \code{FALSE}.
#' @return Outputs a list with limma eBayes fitting, QC plots and a volcano distribution
#' @importFrom limma voom voomWithQualityWeights plotMA plotSA lmFit eBayes topTable plotMDS
#' @importFrom edgeR DGEList calcNormFactors
#' @examples
#' \dontrun{
#' fitRNA <- mirFit(simDfm, wgt = TRUE)
#' }
#' @export
mirFit <- function(dfm, anno, fileName, count_threshold = "none", wgt = FALSE){
  # load annoation
  SampleIndex <- read.csv(file = anno, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  names(SampleIndex)[2] <- "Condition"

  Nrm <- DGEList(counts = dfm[, -1], genes = dfm[1])
  # filtering
  if (count_threshold != "none"){ # set the count threshold for filtering
    count_s <- rowSums(Nrm$counts) # thresholdd
    isexpr <- count_s > count_threshold

    Nrm <- Nrm[isexpr, , keep.lib.size = FALSE] # filtering
  }

  # normalization
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
  with(exportTable, plot(logFC, -log10(P.Value), pch = 20, cex = 0.5, main = "Volcano plot"))

  # Add colored points: red if padj < 0.05, orange of log2FC > 1, green if both)
  with(subset(exportTable, adj.P.Val < 0.05 ), points(logFC, -log10(P.Value), pch = 20, cex = 0.5, col = "red"))
  with(subset(exportTable, abs(logFC) >= log2(1.5)), points(logFC, -log10(P.Value), pch = 20, cex = 0.5,
                                                            col = "orange"))
  with(subset(exportTable, adj.P.Val < 0.05 & abs(logFC) >= log2(1.5)), points(logFC, -log10(P.Value),
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
