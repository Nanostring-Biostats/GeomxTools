setValidity2("NanoStringGeomxSet",
function(object)
{
  msg <- NULL
  if (dim(object)[["Samples"]] > 0L) {
    # sampleNames
    if (!all(grepl("\\.dcc$", sampleNames(object), ignore.case = TRUE))) {
      msg <- c(msg, "'sampleNames' must all have an \".DCC\" file extension")
    }
    # protocolData
    protocolDataColNames <- rownames(.dccMetadata[["protocolData"]])
    if (!all(protocolDataColNames %in% varLabels(protocolData(object)))) {
      msg <-
        c(msg,
          sprintf("'protocolData' must contain columns %s",
                  paste0("\"", protocolDataColNames, "\"", collapse = ", ")))
    }
    # protocolData - FileVersion
    if (!all(protocolData(object)[["FileVersion"]] %in%
             numeric_version(c("0.1")))) {
      msg <-
        c(msg, "'protocolData' \"FileVersion\" must all be 0.1")
    }
  }
  if (dim(object)[["Features"]] > 0L) {
    # featureData
    featureDataColNames <- c("Target")
    if (!all(featureDataColNames %in% varLabels(featureData(object)))) {
      msg <-
        c(msg,
          sprintf("'featureData' must contain columns %s",
                  paste0("\"", featureDataColNames, "\"", collapse = ", ")))
    } 
  }
  if (sum(dim(object)) > 0L) {
    # annotation
    if (length(annotation(object)) != 1L || is.na(annotation(object)) ||
        !nzchar(annotation(object))) {
      msg <- c(msg, "'annotation' must contain the PKC")
    }
  }
  if (prod(dim(object)) > 0L) {
    # assayData
    if (!.validPositiveNumber(exprs(object))) {
      msg <- c(msg, "'exprs' does not contain positive values")
    }
    # dimLabels
    if (length(dimLabels(object)) != 2L) {
      msg <- c(msg, "dimLabels must be a character vector of length 2")
    } else {
      if (!(dimLabels(object)[1L] %in% fvarLabels(object))) {
        msg <- c(msg, "dimLabels[1] must be in 'fvarLabels'")
      }
      if (!(dimLabels(object)[2L] %in% svarLabels(object))) {
        msg <- c(msg, "dimLabels[2] must be in 'svarLabels'")
      }
    }
  }
  if (any(duplicated(c(fvarLabels(object), svarLabels(object),
                       assayDataElementNames(object),
                       "signatures", "design")))) {
    msg <-
      c(msg,
        "'fvarLabels', 'svarLabels', 'assayDataElementNames', \"signatures\", and \"design\" must be unique")
  }
  if (length(signatures(object)) > 0L) {
    numTargets <- lengths(signatures(object))
    if (is.null(names(numTargets)) || any(nchar(names(numTargets)) == 0L)) {
      msg <- c(msg, "'signatures' must be a named NumericList")
    }
    if (any(numTargets == 0L)) {
      msg <- c(msg, "'signatures' vectors must be non-empty")
    } else {
      targets <- names(unlist(unname(stats::weights(signatures(object)))))
      if (is.null(targets) || any(nchar(targets) == 0L)) {
        msg <- c(msg, "'signatures' vectors must be named")
      } else if(!all(unique(targets) %in%
                     c("(Intercept)", featureData(object)[["TargetName"]]))) {
        msg <-
          c(msg,
            "'signatures' vectors must be named with values from 'featureData' \"TargetName\"")
      }
    }
  }
  if (!is.null(design(object))) {
    if (length(design(object)) != 2L) {
      msg <- c(msg, "'design' must be NULL or a one-sided formula")
    }
    if (!all(all.vars(design(object)) %in% varLabels(object))) {
      msg <- c(msg, "'design' must reference columns from 'phenoData'")
    }
  }
  if (is.null(msg)) TRUE else msg
})
