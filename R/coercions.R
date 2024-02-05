
#' @name as.Seurat
NULL

#' Convert GeoMxSet Object to SeuratObject
#' 
#' @inheritParams SeuratObject::as.Seurat 
#' @param ident column in GeoMxSet segmentProperties to set as Seurat object's identity class
#' @param normData assay containing normalized data
#' @param coordinates X and Y coordinates of each ROI, format: c(X,Y)
#' @param forceRaw should raw data be forced into SeuratObject, not recommended
#' @param ... Arguments passed to other methods
#' 
#' @return SeuratObject containing GeoMx data
#' 
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data", package = "GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' 
#' target_demoData <- aggregateCounts(demoData[1:1000,1:10])
#' 
#' target_demoData <- normalize(target_demoData, "quant")
#' 
#' seurat_demoData <- as.Seurat(target_demoData, ident = "cell_line",
#'                              normData = "exprs_norm", forceRaw = FALSE)
#' 
#' @importFrom Seurat CreateSeuratObject 
#' @importFrom Seurat AddMetaData 
#' @importFrom data.table as.data.table
#' 
#' @export 
#' @rdname as.Seurat
#' 
#' @method as.Seurat NanoStringGeoMxSet
#' 
#' @seealso \code{\link[SeuratObject:as.Seurat]{SeuratObject::as.Seurat}}
#' 
#' @concept objects

as.Seurat.NanoStringGeoMxSet <- function(x, ident = NULL, normData = NULL, 
                                         coordinates = NULL, 
                                         forceRaw = FALSE, ...){
    
    if (!try(requireNamespace("Seurat", quietly = TRUE))) {
        stop("Please install Seurat from CRAN before converting to a Seurat object")
    }else{
        requireNamespace("Seurat", quietly = TRUE)
    }
    
    
    if(featureType(x) == "Probe"){
        stop("Data must be on Target level before converting to a Seurat Object")
    }
    
    if(is.null(normData)){
        stop("normData must not be NULL, please provide name of normalized counts matrix")
    }
    
    if(!normData %in% assayDataElementNames(x)){
        stop(paste0("The normData name \"", normData, "\" is not a valid assay name. Valid names are: ", 
                    paste(assayDataElementNames(x), collapse = ", ")))
    }
    
    normFactor_names <- "normFactors|qFactors|negFactors|hkFactors|hknormFactors"
    
    if(length(grep(pattern = normFactor_names, names(sData(x)))) == 0 & 
       forceRaw == FALSE){
        stop("It is NOT recommended to use Seurat's normalization for GeoMx data. 
             Normalize using GeomxTools::normalize() or set forceRaw to TRUE if you want to continue with Raw data")
    }else if(length(grep(pattern = normFactor_names, names(sData(x)))) == 0){
      message("Coercing raw data, it is NOT recommended to use Seurat's normalization for GeoMx data.")
    }
    
    
    sequencingMetrics <- c("FileVersion", "SoftwareVersion", "Date", "Plate_ID", 
                           "Well", "SeqSetId", "Raw", "Trimmed", "Stitched", 
                           "Aligned", "umiQ30", "rtsQ30", "DeduplicatedReads", 
                           "NTC_ID", "NTC", "Trimmed (%)", "Stitched (%)", 
                           "Aligned (%)", "Saturated (%)")
    
    QCMetrics <- "QCFlags"
    
    meta <- as.data.frame(as.data.table(sData(x)[,!colnames(sData(x)) %in%  c(sequencingMetrics, QCMetrics)]))
    
    if(any(grepl("_", rownames(x)))){
      rownames(x) <- gsub("_", "-", rownames(x))
      message("Feature names cannot have underscores ('_'), replacing with dashes ('-')")
    }
    
    if(packageVersion("Seurat") < 5){
      seuratConvert <- suppressWarnings(Seurat::CreateSeuratObject(counts = assayDataElement(x, normData), 
                                                                   assay = "GeoMx", 
                                                                   project = expinfo(experimentData(x))[["title"]]))
      seuratConvert <- suppressWarnings(Seurat::AddMetaData(object = seuratConvert, 
                                                            metadata = meta))
      seuratConvert@assays$GeoMx <- Seurat::AddMetaData(object = seuratConvert@assays$GeoMx, 
                                                        metadata = fData(x))
      
      if(!is.null(ident)){
        if(!ident %in% colnames(seuratConvert@meta.data)){
          stop(paste0("ident \"", ident, "\" not found in GeoMxSet Object"))
        }
        
        Seurat::Idents(seuratConvert) <- seuratConvert[[ident]]
      }
    }else{
      projectName <- expinfo(experimentData(x))[["title"]]
      if(projectName == ""){
        projectName <- "GeoMx"
      }
      
      seuratConvert <- suppressWarnings(Seurat::CreateSeuratObject(counts = assayDataElement(x, "exprs"), 
                                                                   assay = "GeoMx", 
                                                                   project = projectName))
      seuratConvert <- Seurat::SetAssayData(seuratConvert, layer = "data", 
                                            new.data = assayDataElement(x, normData))
      seuratConvert <- suppressWarnings(Seurat::AddMetaData(object = seuratConvert, 
                                                            metadata = meta))
      seuratConvert@assays$GeoMx <- Seurat::AddMetaData(object = seuratConvert@assays$GeoMx, 
                                                        metadata = fData(x))
      
      if(!is.null(ident)){
        if(!ident %in% colnames(seuratConvert@meta.data)){
          stop(paste0("ident \"", ident, "\" not found in GeoMxSet Object"))
        }
        
        Seurat::Idents(seuratConvert) <- as.factor(seuratConvert@meta.data[[ident]])
      }
    }

    seuratConvert@misc <- otherInfo(experimentData(x)) 
    seuratConvert@misc[["sequencingMetrics"]] <- sData(x)[colnames(sData(x)) %in% 
                                                                   sequencingMetrics]
    seuratConvert@misc[["QCMetrics"]] <- sData(x)[colnames(sData(x)) %in% 
                                                           QCMetrics]
    
    if(ncol(seuratConvert@misc[["QCMetrics"]]) == 0){
      seuratConvert@misc[["QCMetrics"]] <- NULL
    }
    
    if(!is.null(coordinates)){
        xcoord <- coordinates[1]
        ycoord <- coordinates[2]
        
        if(xcoord %in% colnames(seuratConvert@meta.data) & 
           ycoord %in% colnames(seuratConvert@meta.data)){
            coord.df <- data.frame(x=seuratConvert@meta.data[[xcoord]], 
                                   y=seuratConvert@meta.data[[ycoord]])
            colnames(coord.df) <- coordinates
            seuratConvert@meta.data <- seuratConvert@meta.data[!colnames(seuratConvert@meta.data) %in% 
                                                                   coordinates]
        }else{
            if(!xcoord %in% colnames(seuratConvert@meta.data) &
               !ycoord %in% colnames(seuratConvert@meta.data)){
                stop(paste0("xcoord \"", xcoord, "\" and ycoord \"", 
                            ycoord, "\" not found in GeoMxSet Object"))
            }
            
            if(!xcoord %in% colnames(seuratConvert@meta.data)){
                stop(paste0("xcoord \"", xcoord, 
                            "\" not found in GeoMxSet Object"))
            }
   
            if(!ycoord %in% colnames(seuratConvert@meta.data)){
                stop(paste0("ycoord \"", ycoord, 
                            "\" not found in GeoMxSet Object"))
            }
        }
        
        rownames(coord.df) <- rownames(seuratConvert@meta.data)
        
        # need to create DSP specific image class
        seuratConvert@images$image =  new(
            Class = 'SlideSeq',
            assay = "GeoMx",
            key = "image_",
            coordinates = coord.df
        )
    }
    
    return(seuratConvert)
}

#' Convert Object to SpatialExperiment
#' 
#' @param x GeoMxSet object to convert 
#' @param ... arguments to be passed to other methods
#' 
#' @export

as.SpatialExperiment <- function(x, ...) {
    UseMethod(generic = 'as.SpatialExperiment', object = x)
}

#' Convert GeoMxSet Object to SpatialExperiment
#' 
#' @param x GeoMxSet object to convert 
#' @param normData assay containing normalized data
#' @param coordinates X and Y coordinates of each ROI, format: c(X,Y)
#' @param forceRaw should raw data be forced into SpatialExperiment, not recommended
#' @param ... Arguments passed to other methods
#' 
#' @return SpatialExperiment containing GeoMx data
#' 
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data", 
#'                        package = "GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' 
#' target_demoData <- aggregateCounts(demoData[1:1000,1:10])
#' 
#' target_demoData <- normalize(target_demoData, "quant")
#' 
#' seurat_demoData <- as.SpatialExperiment(target_demoData, 
#'                                         normData = "exprs_norm", 
#'                                         forceRaw = FALSE)
#' 
#' @importFrom SpatialExperiment SpatialExperiment
#' 
#' @export 
#' 
#' @rdname as.SpatialExperiment
#' @method as.SpatialExperiment NanoStringGeoMxSet

as.SpatialExperiment.NanoStringGeoMxSet <- function(x, normData = NULL, 
                                                    coordinates = NULL, 
                                                    forceRaw = FALSE, ...){
    
    if (!try(requireNamespace("SpatialExperiment", quietly = TRUE))) {
        stop("Please install SpatialExperiment from Bioconductor before converting to a SpatialExperiment object")
    }else{
        requireNamespace("SpatialExperiment", quietly = TRUE)
    }
    
    if(featureType(x) == "Probe"){
        stop("Data must be on Target level before converting to a SpatialExperiment Object")
    }
    
    if(is.null(normData)){
        stop("normData must not be NULL, please provide name of normalized counts matrix")
    }
    
    if(!normData %in% assayDataElementNames(x)){
        stop(paste0("The normData name \"", normData, 
                    "\" is not a valid assay name. Valid names are: ", 
                    paste(assayDataElementNames(x), collapse = ", ")))
    }
    
    normFactor_names <- "normFactors|qFactors|negFactors|hkFactors|hknormFactors"
    
    if(length(grep(pattern = normFactor_names, names(sData(x)))) == 0 & 
       forceRaw == FALSE){
        stop("It is NOT recommended to use Seurat's normalization for GeoMx data. 
             Normalize using GeomxTools::normalize() or set forceRaw to TRUE if you want to continue with Raw data")
    }else if(length(grep(pattern = normFactor_names, names(sData(x)))) == 0){
      warning("Coercing raw data, it is NOT recommended to use Seurat's normalization for GeoMx data.")
    }
    
    sequencingMetrics <- c("FileVersion", "SoftwareVersion", "Date", "Plate_ID", 
                           "Well", "SeqSetId", "Raw", "Trimmed", "Stitched", 
                           "Aligned", "umiQ30", "rtsQ30", "DeduplicatedReads", 
                           "NTC_ID", "NTC", "Trimmed (%)", "Stitched (%)", 
                           "Aligned (%)", "Saturated (%)")
    
    QCMetrics <- "QCFlags"
    
    if(!is.null(coordinates)){
        xcoord <- coordinates[1]
        ycoord <- coordinates[2]
        
        if(xcoord %in% colnames(sData(x)) & 
           ycoord %in% colnames(sData(x))){
            coord.df <- data.frame(x=sData(x)[[xcoord]], 
                                   y=sData(x)[[ycoord]])
            colnames(coord.df) <- coordinates
            sequencingMetrics <- c(sequencingMetrics, coordinates)
        }else{
            if(!xcoord %in% colnames(sData(x)) &
               !ycoord %in% colnames(sData(x))){
                stop(paste0("xcoord \"", xcoord, "\" and ycoord \"", ycoord, 
                            "\" not found in GeoMxSet Object"))
            }
            
            if(!xcoord %in% colnames(sData(x))){
                stop(paste0("xcoord \"", xcoord, 
                            "\" not found in GeoMxSet Object"))
            }
            
            if(!ycoord %in% colnames(sData(x))){
                stop(paste0("ycoord \"", ycoord, 
                            "\" not found in GeoMxSet Object"))
            }
        }
        
        rownames(coord.df) <- rownames(sData(x))
        
        coord.df <- as.matrix(coord.df)
        
    }else{
        coord.df <- NULL
    }
    
    spe <- SpatialExperiment::SpatialExperiment(assay = assayDataElement(x, normData),
                                                colData = sData(x)[!colnames(sData(x)) %in% 
                                                                            c(sequencingMetrics, 
                                                                              QCMetrics)],
                                                rowData = fData(x),
                                                spatialCoords = coord.df)
    
    names(spe@assays@data) <- "GeoMx"
    
    spe@metadata <- otherInfo(experimentData(x))
    spe@metadata[["sequencingMetrics"]] <- sData(x)[colnames(sData(x)) %in% 
                                                             sequencingMetrics]
    spe@metadata[["QCMetrics"]] <- sData(x)[colnames(sData(x)) %in% 
                                                     QCMetrics]
    
    if(ncol(spe@metadata[["QCMetrics"]]) == 0){
      spe@metadata[["QCMetrics"]] <- NULL
    }
    
    return(spe)
}

#' Update GeoMxSet object to current version
#' 
#' @param object GeoMxSet object to update 
#' 
#' @return updated GeoMxSet object
#' 
#' @importFrom BiocGenerics getObjectSlots
#' 
#' @export 
#' 
updateGeoMxSet <- function(object){
  if(!"analyte" %in% names(getObjectSlots(object))){
    object@analyte <- "RNA"
    
    object@.__classVersion__$NanoStringGeoMxSet <- 
        paste(packageVersion("GeomxTools"), collapse=".")
  }else{
    warning("GeoMxSet object up to date, no update necessary")
  }

  return(object)
}

