#' Produce concordance pairs plots
#' 
#' @description 
#' Draws a pairs plot of concordance between multiple factors/ probes
#' Upper panels are the concordance plot.
#' Lower panels are the standard deviation of the log2-ratios between the targets 
#' 
#' @param dat dataframe of values to be compared. Could be e.g. IgG probe counts,
#'  or different normalization factors. Each variable/factor is in a column.
#' @param plotFactor segment factor to color the plot by
#' 
#' @importFrom expr_text rlang
#' @importFrom ggpairs GGally
#' @importFrom wrap GGally
#' 
#' @noRd
#' 
concordancePlot <- function(dat, plotFactor){
    n <- ncol(dat)
    p <- GGally::ggpairs(dat, mapping = ggplot2::aes(colour=dat[[plotFactor]], 
                                                     alpha=0.1),
                         columns=1:n, progress=FALSE,
                         lower = list(continuous = GGally::wrap("smooth",
                                                                alpha = 0.3,
                                                                size=0.1)),
                         upper = list(continuous = print.sd.log.ratio)) +
        ggplot2::theme_bw()
    
    
    return(p)
}

#' Prints the standard deviation of log2-ratios
#' 
#' @description 
#' Similar calculation to correlation but is not affected by high differences 
#' across values. 
#' 
#' @param data dataframe of XY values
#' @param mapping ggplot mapping value
#' 
#' @noRd
print.sd.log.ratio <- function(data, mapping) {
  x <- gsub("\\n", "\n", 
            gsub("~", "", 
                 gsub(pattern = "`", "", 
                      rlang::expr_text(mapping$x))), 
            fixed = TRUE)
  y <- gsub("\\n", "\n", 
            gsub("~", "", 
                 gsub(pattern = "`", "", 
                      rlang::expr_text(mapping$y))), 
            fixed = TRUE)
  
  x <- data[[x]]
  y <- data[[y]]
  
  center <- function(x){(max(x) + min(x))/2}
  
  pairwise.stat <- sd(log2(pmax(x, 1)) - log2(pmax(y, 1)), na.rm = TRUE)
  
  return(ggplot2::ggplot(data) +
           ggplot2::geom_blank(mapping) +
           ggplot2::annotate(geom = "text", x = center(x), y = center(y),
                             label = paste("SD(log2(ratios)) =\n", 
                                           round(pairwise.stat, 3))))
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
    ylim <- range(log2(snr))
    if(ylim[1L] == -Inf){
      ylim[1L] <- 0
    }
    if(ylim[2L] == -Inf){
      ylim[2L] <- 0
    }
    
    fig <- function(){
        par(mar = c(11, 4, 2, 1))
        boxplot(t(log2(snr)),
                las = 2,
                outline = FALSE,
                ylim = ylim,
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
    }
    
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
#' @noRd

snrOrder <- function(object, neg.names){
    if(analyte(object) != "Protein"){
        stop("This function is only meant for protein data")
    }  
    
    if(is.null(neg.names)){
        neg.names <- iggNames(object)
    }
    
    raw <- exprs(object)
    
    # estimate background:
    negfactor <- apply(raw[neg.names, , drop = FALSE], 2, 
                       function(x){pmax(mean(x), 1)})
    
    # calc snr
    snr <- sweep(raw, 2, negfactor, "/")
    
    igginds <- which(is.element(rownames(snr), neg.names))
    o <- c(igginds, setdiff(order(apply(snr, 1, median)), igginds))
    
    return(snr[o,, drop=FALSE])
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
    return(rownames(snrOrder(object, neg.names)))
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
#' protSegTypeFig
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
    
    
    tempmat <- t(pmax(exprs(object)[targetList, ], 1))
    colnames(tempmat) <- paste0(colnames(tempmat), " counts")
    
    fig <- concordancePlot(dat = as.data.frame(cbind(tempmat, 
                                                     sData(object)[as.character({{plotFactor}})])),
                           plotFactor = plotFactor)
    
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
#' normConcord
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
    
    if(is.null(normfactors)){
        normfactors <- computeNormalizationFactors(object)
    }
    
    if(ncol(normfactors) <= 1){
        stop("At least 2 normfactors must be given for comparison")
    }
    
    # pairs plots:
    tempmat <- pmax(normfactors, 1)
    colnames(tempmat)[colnames(tempmat) == "HK geomean"] <- "HK geomean\n(counts)"
    colnames(tempmat)[colnames(tempmat) == "Neg geomean"] <- "Neg geomean\n(counts)"
    colnames(tempmat)[colnames(tempmat) == "Nuclei"] <- "Nuclei\n(counts)"
    colnames(tempmat)[colnames(tempmat) == "Area"] <- "Area\n(microns squared)"
    
    fig <- concordancePlot(dat = as.data.frame(cbind(tempmat, 
                                                     sData(object)[as.character({{plotFactor}})])),
                           plotFactor = plotFactor)
    
    return(fig)
}

