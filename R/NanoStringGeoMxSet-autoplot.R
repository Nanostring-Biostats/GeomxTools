#' Produce concordance pairs plots
#' 
#' @description 
#' Draws a pairs plot of concordance between multiple factors/ probes
#' Upper panels are the concordance plot.
#' Lower panels are the standard deviation of the log2-ratios between the targets 
#' 
#' @param mat Matrix of values to be compared. Could be e.g. IgG probe counts,
#'  or different normalization factors. Each variable/factor is in a column.
#' @param col Vector of colors for points
#' @param collegend Named vector of colors, used to draw a legend.
#' @param legend.main Title for the color legend, used to name the variable of interest
#' @param main Plot title
#' @param pch point type argument passed to pairs()
#' @param cex point size argument passed to pairs()
#' @param ... Additional arguments passed to pairs()
#' @examples
#' # simulate HKs:
#' x = pmax(rnorm(100, 50, 10), 0)
#' mat = sweep(matrix(rnorm(300, 0, 5), 100), 1, x, "+")
#' concordancePlot(mat = mat, col = rep(c("blue", "orange"), each = 50))

concordancePlot <- function(mat, col = rgb(0, 0, 0, 0.5),
                             collegend = NULL, legend.main = NULL,
                             main = "", pch = 16, cex = 1.5, ...) {
  
  # subsidiary function for printing the SD of a ratio, for use by pairs():
  print.sd.log.ratio <- function(x, y, digits = 2, prefix = "", cex.cor = 1, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    pairwise.stat <- sd(log2(pmax(x, 1)) - log2(pmax(y, 1)), na.rm = TRUE)
    txt <- round(pairwise.stat, 2)
    txt <- paste0(prefix, txt)
    legend("center", legend = txt, box.col = "white", cex = 1.5)
  }
  
  # draw pairs plot:
  par(mar = c(4, 4, 2, 1))
  
  concordance_plot <- function(){
    pairs(mat,
          log = "xy",
          upper.panel = points,
          lower.panel = print.sd.log.ratio,
          oma = c(3, 3, 3, 35),
          col = col,
          pch = pch,
          cex = cex,
          labels = colnames(mat),
          main = main
    )
    # draw a legend:
    if (length(collegend) > 0) {
      legend("right",
             col = c(NA, collegend, NA, NA),
             legend = c(legend.main, names(collegend), "", "Numbers show SD(log2(ratios))"),
             pch = 16
      )
    }
  }
  
  return(concordance_plot)
}

#' Generate Protein QC signal boxplot figure
#' 
#' @description 
#' For use with protein data ONLY.
#' 
#' @param object name of the object class to subset
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
#' }
#' @param neg.names names of IgGs, if NULL IgGs will be detected automatically
#' 
#' @return figure function
#' 
#' @examples
#' proteinData <- readRDS(file= system.file("extdata","DSP_Proteogenomics_Example_Data", 
#' "proteinData.rds", package = "GeomxTools"))
#' 
#' igg.names <- iggNames(proteinData)
#' 
#' qcFig <- qcProteinSignal(object = proteinData, neg.names = igg.names)
#' 
#' qcFig()
#' 
#' @export

qcProteinSignal <- function(object, neg.names=NULL) {
  if(analyte(object) != "Protein"){
    stop("This figure is only meant for protein data")
  }  
  
  if(is.null(neg.names)){
    neg.names <- iggNames(object)
  }
  
  snr <- snrOrder(object, neg.names)
  
  protnames <- rownames(snr)
  
  fig <- function(){
    par(mar = c(11, 4, 2, 1))
    boxplot(t(log2(snr)),
            las = 2,
            outline = FALSE,
            ylim = range(log2(snr)),
            names = protnames,
            ylab = "Log2 signal-to-background ratio",
            cex.axis = .85 - 0.3 * (nrow(snr) > 60)
    )
    axis(2, at = 1, labels = 1, las = 2, cex = 0.5)
    points(jitter(rep(1:nrow(snr), ncol(snr))),
           log2(snr),
           col = "#00008B80",
           pch = 16, cex = 0.5
    )
    abline(h = 0)
    abline(v = length(neg.names) + 0.5, lty = 2)
    abline(h = 1, lty = 2)
  }
  
  fig()
  
  return(fig)
}

#' Generate ordered SNR matrix
#' 
#' @description 
#' For use with protein data ONLY.
#' 
#' @param object name of the object class to subset
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
#' }
#' @param neg.names names of IgGs, if NULL IgGs will be detected automatically
#' 
#' @return SNR matrix in increasing order
#' 

snrOrder <- function(object, neg.names){
  if(analyte(object) != "Protein"){
    stop("This figure is only meant for protein data")
  }  
  
  if(is.null(neg.names)){
    neg.names <- iggNames(object)
  }
  
  raw <- exprs(object)
  
  # estimate background:
  negfactor <- apply(raw[neg.names, , drop = FALSE], 2, function(x){pmax(ngeoMean(x), 1)})
  
  # calc snr
  snr <- sweep(raw, 2, negfactor, "/")
  
  igginds <- which(is.element(rownames(snr), neg.names))
  o <- c(igginds, setdiff(order(apply(snr, 1, median)), igginds))
  
  return(snr[o,])
}

#' Generate list of proteins ordered by SNR
#' 
#' @description 
#' For use with protein data ONLY.
#' 
#' @param object name of the object class to subset
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
#' }
#' @param neg.names names of IgGs, if NULL IgGs will be detected automatically
#' 
#' @return protein names in increasing SNR order
#' 
#' @examples
#' proteinData <- readRDS(file= system.file("extdata","DSP_Proteogenomics_Example_Data", 
#' "proteinData.rds", package = "GeomxTools"))
#' 
#' igg.names <- iggNames(proteinData)
#' 
#' proteinOrder <- qcProteinSignalNames(object = proteinData, neg.names = igg.names)
#' 
#' @export

qcProteinSignalNames <- function(object, neg.names){
  if(analyte(object) != "Protein"){
    stop("This figure is only meant for protein data")
  }  
  
  if(is.null(neg.names)){
    neg.names <- iggNames(object)
  }
  
  raw <- exprs(object)
  
  # estimate background:
  negfactor <- apply(raw[neg.names, , drop = FALSE], 2, function(x){pmax(ngeoMean(x), 1)})
  
  # calc snr
  snr <- sweep(raw, 2, negfactor, "/")
  
  igginds <- which(is.element(rownames(snr), neg.names))
  o <- c(igginds, setdiff(order(apply(snr, 1, median)), igginds))
  
  return(rownames(snr[o,]))
}

#' Generate concordance figure of targets based on user provided factors
#' 
#' @description  
#' Upper panels are the concordance plot.
#' Lower panels are the standard deviation of the log2-ratios between the targets 
#' 
#' @param targetList names of targets to plot concordance, normally IgGs. 
#' @param object name of the object class to subset
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
#' }
#' @param plotFactor segment factor to color the plot by
#' 
#' @examples
#' proteinData <- readRDS(file= system.file("extdata","DSP_Proteogenomics_Example_Data", 
#' "proteinData.rds", package = "GeomxTools"))
#' 
#' igg.names <- iggNames(proteinData)
#' 
#' protSegTypeFig <- plotConcordance(targetList = igg.names, object = proteinData, 
#'                                   plotFactor = "Segment_Type")
#' protSegTypeFig()
#' 
#' RNASegTypeFig <- plotConcordance(targetList = c("C1orf43", "GPI", "OAZ1"), 
#'                                  object = RNAData,plotFactor = "Segment_Type")
#' RNASegTypeFig()
#' 
#' @export

plotConcordance <- function(targetList, object, plotFactor){
  
  if(!plotFactor %in% colnames(sData(object))){
    stop("Given plotFactor are not in dataset, spelling and capitalization matter")
  }
  
  if (length(targetList) <= 1) {
    stop("At least 2 targets must be given for comparisons")
  }
  
  if(length(plotFactor) > 1){
    plotFactor <- plotFactor[1]
    warning("Only the first plotFactor will be plotted, please call function again for other factors")
  }
  
  if(all(!targetList %in% rownames(object))){
    stop("Given targetList does not match target names in object")
  }else if(any(!targetList %in% rownames(object))){
    notIn <- targetList[which(!targetList %in% rownames(object))]
    warning(paste("These targets from targetList do not match target names in object and will not be part of analysis:", 
                  paste(notIn, collapse = ", ")))
  }
  
  cols <- assignColors(annot = sData(object)[, plotFactor, drop = FALSE])
  
  par(mar = c(4, 4, 4, 1))

  tempmat <- t(pmax(exprs(object)[targetList, ], 1))
  colnames(tempmat) <- paste0(colnames(tempmat), " counts")
  fig <- concordancePlot(
    mat = tempmat,
    col = cols[[plotFactor]][as.character(sData(object)[, plotFactor])],
    collegend = cols[[plotFactor]],
    legend.main = plotFactor
  )
  
  fig()

  return(fig)
}


#' Generate concordance figure of normalization factors based on user provided factors
#' 
#' @description 
#' For use with protein data ONLY.
#' 
#' Upper panels are the concordance plot.
#' Lower panels are the standard deviation of the log2-ratios between the normalization factors 
#' 
#' @param object name of the object class to subset
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
#' }
#' @param plotFactor segment factor to color the plot by
#' @param normfactors normalization factors from computeNormalizationFactors(). If NULL these are calculated automatically. 
#' 
#' @examples
#' proteinData <- readRDS(file= system.file("extdata","DSP_Proteogenomics_Example_Data", 
#' "proteinData.rds", package = "GeomxTools"))
#' 
#' normConcord <- plotNormFactorConcordance(object = proteinData, plotFactor = "Segment_Type")
#' normConcord()
#' 
#' @export

plotNormFactorConcordance <- function(object, plotFactor, normfactors = NULL){
  
  if(analyte(object) != "Protein"){
    stop("This figure is only meant for protein data")
  }  
  
  if(length(plotFactor) > 1){
    plotFactor <- plotFactor[1]
    warning("Only the first plotFactor will be plotted, please call function again for other factors")
  }
  
  if(!plotFactor %in% colnames(sData(object))){
    stop("Given plotFactor are not in dataset, spelling and capitalization matter")
  }
  
  cols <- assignColors(annot = sData(object)[, plotFactor, drop = FALSE])
  
  if(is.null(normfactors)){
    normfactors <- computeNormalizationFactors(object)
  }
  
  if(ncol(normfactors) <= 1){
    stop("At least 2 normfactors must be given for comparison")
  }
  
  par(mar = c(4, 4, 2, 1))
  # pairs plots:
  tempmat <- pmax(normfactors, 1)
  colnames(tempmat)[colnames(tempmat) == "HK geomean"] <- "HK geomean\n(counts)"
  colnames(tempmat)[colnames(tempmat) == "Neg geomean"] <- "Neg geomean\n(counts)"
  colnames(tempmat)[colnames(tempmat) == "Nuclei"] <- "Nuclei\n(counts)"
  colnames(tempmat)[colnames(tempmat) == "Area"] <- "Area (microns squared)"
  
  fig <- concordancePlot(
    mat = tempmat,
    col = cols[[plotFactor]][as.character(sData(object)[, plotFactor])],
    collegend = cols[[plotFactor]],
    legend.main = plotFactor,
    main = "Normalization factors"
  ) 
  
  fig()
  
  return(fig)
}

assignColors <- function(annot) {
  # vector of colors to choose from:
  colvec <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
    "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", 
    "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", 
    "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", 
    "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#A6CEE3", 
    "#B2DF8A", "#FB9A99", "#FDBF6F", "#CAB2D6", "#FFFF99", "#33A02C", 
    "#FF7F00", "#B15928", "#1F78B4", "#E31A1C", "#6A3D9A", sample(colors(), 100)
  )
  # subsidiary function to fade colors (works like the "alpha" function from the "scales" package)
  fadecols <- function(cols, fade = .5) {
    fcols <- c()
    for (i in 1:length(cols))
    {
      tmp <- as.vector(col2rgb(cols[i]) / 256)
      fcols[i] <- rgb(tmp[1], tmp[2], tmp[3], fade)
    }
    return(fcols)
  }
  colvec <- fadecols(colvec, 0.7)
  
  # assign colors:
  cols <- list()
  colorby <- colnames(annot)
  for (varname in colorby) {
    varlevels <- as.character(unique(annot[, varname]))
    cols[[varname]] <- colvec[1:length(varlevels)]
    names(cols[[varname]]) <- varlevels
    # remove the used colors from further consideration:
    colvec <- setdiff(colvec, cols[[varname]]) # (disabling this so the more bold colors are re-used)
  }
  
  return(cols)
}
