#' Run a mixed model on GeoMxSet
#'
#' @param object name of the object class to perform QC on
#' \enumerate{
#'     \item{NanoStringGeoMxSet, use the NanoStringGeoMxSet class}
#' }
#' @param elt assayDataElement of the geoMxSet object to run the DE on
#' @param modelFormula formula used in DE, if null, the design(object) is used
#' @param groupVar = "group",  sample annotation to group the data for comparing means
#' @param nCores = 1, number of cores to use, set to 1 if running in serial mode
#' @param multiCore = TRUE, set to TRUE to use multiCore, FALSE to run in cluster mode
#' @param pAdjust = "BY" method for p-value adjustment
#' @param pairwise boolean to calculate least-square means pairwise differences
#'
#' @return mixed model output list
#'
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data", package = "GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' target_demoData <- aggregateCounts(demoData)
#' target_demoData <- normalize(target_demoData, norm_method="quant")
#' target_demoData <- target_demoData[1:100, ]
#' pData(target_demoData)[["slide"]] <- 
#'     factor(pData(target_demoData)[["slide name"]])
#' protocolData(target_demoData)[["pool_rep"]] <- 
#'     factor(protocolData(target_demoData)[["pool_rep"]])
#' mixedOutmc <- mixedModelDE(target_demoData,
#'                            elt = "exprs_norm",
#'                            modelFormula = ~ pool_rep +  (1 | slide),
#'                            groupVar = "pool_rep",
#'                            nCores = 12,
#'                            multiCore = TRUE,
#'                            pAdjust = NULL
#' )
#'
#' @export
#'

mixedModelDE <- function(object, elt = "exprs", modelFormula = NULL,
                         groupVar = "group", nCores = 1, multiCore = TRUE,
                         pAdjust = "BY", pairwise = TRUE) {
  if (is.null(modelFormula)) {
    modelFormula <- design(object)
  }
  mTerms <- all.vars(modelFormula)
  if ("1" %in% mTerms) {
    mTerms <- mTerms[which(!(mTerms %in% "1"))]
  }
  # check if groupVar is in model formula terms
  if (!groupVar %in% mTerms){
    stop ("Error: groupVar needs to be defined as fixed effect in the model.\n")
  }
  # check if terms in model are in sData
  if (any(!mTerms %in% names(sData(object)))){
    stop ("Error: Not all terms in the model formula are in pheno or protocol data.\n")
  }
  pDat <- sData(object)[,mTerms]
  for (i in names(pDat))
  {
    if (inherits(i, "character")) {
      pDat[, i] <- as.factor(pDat[, i])
    }
  }
  if (nCores > 1) {
    deFunc <- function(i, groupVar, pDat, modelFormula, exprs, pairwise = TRUE) {
      dat <- data.frame(expr = exprs$exprs[i, ], pDat)
      lmOut <- suppressWarnings(lmerTest::lmer(modelFormula, dat))
      if(pairwise == FALSE) {
        lsm <- lmerTest::ls_means(lmOut, which = groupVar, pairwise = FALSE)
      } else {
        lsm <- lmerTest::ls_means(lmOut, which = groupVar, pairwise = TRUE)
      }
      lmOut <- matrix(stats::anova(lmOut)[groupVar, "Pr(>F)"], ncol = 1, dimnames = list(groupVar, "Pr(>F)"))
      lsmOut <- matrix(cbind(lsm[,"Estimate"], lsm[,"Pr(>|t|)"]), ncol = 2, dimnames = list(gsub(groupVar, "", rownames(lsm)), c("Estimate", "Pr(>|t|)")))

      return(list(anova = lmOut, lsmeans = lsmOut))
    }
    exprs <- new.env()
    exprs$exprs <- assayDataElement(object, elt = elt)
    if (multiCore & Sys.info()['sysname'] != "Windows") {
      mixedOut <- parallel::mclapply(featureNames(object), deFunc, groupVar, pDat, formula(paste("expr", as.character(modelFormula)[2], sep = " ~ ")), exprs, mc.cores = nCores)
    }
    else {
      cl <- parallel::makeCluster(getOption("cl.cores", nCores))
      mixedOut <- parallel::parLapply(cl, featureNames(object), deFunc, groupVar, pDat, formula(paste("expr", as.character(modelFormula)[2], sep = " ~ ")), exprs, pairwise)
      suppressWarnings(parallel::stopCluster(cl))
    }
    mixedOut <- rbind(array(lapply(mixedOut, function(x) x[["anova"]])),
                      array(lapply(mixedOut, function(x) x[["lsmeans"]])))
    colnames(mixedOut) <- featureNames(object)
    rownames(mixedOut) <- c("anova", "lsmeans")
  }
  else {
    deFunc <- function(expr, groupVar, pDat, modelFormula, pairwise = TRUE) {
      dat <- data.frame(expr = expr, pDat)
      lmOut <- suppressMessages(lmerTest::lmer(modelFormula, dat))
      if(pairwise == FALSE) {
        lsm <- lmerTest::ls_means(lmOut, which = groupVar, pairwise = FALSE)
      } else {
        lsm <- lmerTest::ls_means(lmOut, which = groupVar, pairwise = TRUE)
      }
      lmOut <- matrix(stats::anova(lmOut)[groupVar, "Pr(>F)"], ncol = 1, dimnames = list(groupVar, "Pr(>F)"))
      lsmOut <- matrix(cbind(lsm[,"Estimate"], lsm[,"Pr(>|t|)"]), ncol = 2, dimnames = list(gsub(groupVar, "", rownames(lsm)), c("Estimate", "Pr(>|t|)")))

      return(list(anova = lmOut, lsmeans = lsmOut))
    }
    mixedOut <- assayDataApply(object, 1, deFunc, groupVar, pDat, formula(paste("expr", as.character(modelFormula)[2], sep = " ~ ")), pairwise,  elt = elt)
  }
  if (!is.null(pAdjust)) {
    mixedOut["anova", ] <- p.adjust(mixedOut["anova", ], method = pAdjust)
  }
  return(mixedOut)
}
