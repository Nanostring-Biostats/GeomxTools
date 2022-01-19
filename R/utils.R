#' Get the geometric mean of a vector
#' 
#' @param x numeric vector
#' @param thresh minimum numeric value greater than 0 to have in vector
#' 
#' @return numeric geometric mean of vector
#' 
#' @examples
#' ngeoMean(c(0, 1, 2, 2), thresh=0.1)
#' 
#' @export
#' 
ngeoMean <- function(x, thresh=0.5) {
    x <- thresholdValues(x, thresh=thresh)
    return(EnvStats::geoMean(x, na.rm = TRUE))
}

#' Get the geometric standard deviation of a vector
#' 
#' @param x numeric vector
#' @param thresh minimum numeric value greater than 0 to have in vector
#' 
#' @return numeric geometric standard deviation of vector
#' 
#' @examples
#' ngeoSD(c(0, 1, 2, 2), thresh=0.1)
#' 
#' @export
#' 
ngeoSD <- function(x, thresh=0.5) {
    x <- thresholdValues(x, thresh=thresh)
    return(EnvStats::geoSD(x, na.rm=T))
}

#' Get take the log of a numeric vector
#' 
#' @param x numeric vector
#' @param thresh minimum numeric value greater than 0 to have in vector
#' @param base numeric value indicating base to log with
#' 
#' @return numeric vector with logged values
#' 
#' @examples
#' logtBase(c(0, 1, 2, 2), thresh=0.1, base=10)
#' 
#' @export
#' 
logtBase <- function(x, thresh=0.5, base=2) {
    x <- thresholdValues(x, thresh=thresh)
    return(log(x, base=base))
}

thresholdValues <- function(x, thresh=0.5) {
    if (thresh <= 0) {
      warning("Parameter, thresh, cannot be set to less than or equal to 0. 
              The default threshold of 0.5 was used instead.")
      thresh <- 0.5
    }
    if (min(x, na.rm = TRUE) < thresh) {
        x <- x + thresh
    }
    return(x)
}

#' Add one to all counts in an expression matrix
#' 
#' @param object name of the NanoStringGeoMxSet object
#' @param elt expression matrix element in \code{assayDataElement}
#'        to shift all counts by
#' @param useDALogic boolean to use the same logic in DA (impute 0s to 1s)
#'        or set to FALSE to shift all counts by 1
#' 
#' @return object of NanoStringGeoMxSet class
#' 
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data",
#'                        package="GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' shiftCountsOne(demoData)
#' 
#' @export
#' 
shiftCountsOne <- function(object, elt="exprs", useDALogic=FALSE) {
    if (countsShiftedByOne(object)) {
        stop("The exprs matrix has already been shifted by one. ",
             "This operation will not be repeated.")
    }
    assayDataElement(object, "rawZero") <- 
        assayDataElement(object, elt=elt)
    experimentData(object)@other$shiftedByOne <- TRUE
    if (useDALogic) {
        assayDataElement(object, 
                         "exprs")[assayDataElement(object, "exprs") < 1] <- 1
        experimentData(object)@other$shiftedByOneLogic <- 
            "Only zeroes increased by 1 count in exprs matrix"
    } else {
        assayDataElement(object, "exprs") <- 
            assayDataElement(object, elt=elt) + 1
        experimentData(object)@other$shiftedByOneLogic <- 
            "All counts increased by 1 throughout exprs matrix"
    }
    return(object)
}



#' Return the IgG negative controls for protein
#' 
#' @param object name of the NanoStringGeoMxSet object
#' 
#' @return names of IgGs
#' 
#' @export
#' 
igg_names <- function(object){
  names <- featureData(object)$Target[featureData(object)$AnalyteType == "Protein" &
                                        featureData(object)$CodeClass == "Negative"]
  return(names)
}

#' Return the House Keeper positive controls for protein
#' 
#' @param object name of the NanoStringGeoMxSet object
#' 
#' @return names of HKs
#' 
#' @export
#' 
hk_names <- function(object){
  names <- featureData(object)$Target[featureData(object)$AnalyteType == "Protein" &
                                        featureData(object)$CodeClass == "Control"]
  return(names)
}

assign_colors <- function(annot) {
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

#### NOT TESTED OR USED ####
collapseCounts <- function(object) {
    probeCounts <- data.table(cbind(fData(object)[, c("TargetName", "Module")],
                                    assayDataElement(object, elt="exprs")))
    collapsedCounts <- probeCounts[, lapply(.SD, ngeoMean), 
                                     by=c("TargetName", "Module")]
    return(collapsedCounts)
}
