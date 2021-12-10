#' Convert GeoMxSet Object to SeuratObject
#' 
#' @param object GeoMxSet object to convert 
#' @param ident column in GeoMxSet segmentProperties to set as Seurat object's identity class
#' @param normData assay containing normalized data
#' @param coordinates X and Y coordinates of each ROI, format: c(X,Y)
#' @param forceRaw should raw data be forced into SeuratObject, not recommended
#' 
#' @return SeuratObject containing GeoMx data
#' 
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data", package = "GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' 
#' target_demoData <- aggregateCounts(demoData)
#' 
#' target_demoData <- normalize(target_demoData, "quant")
#' 
#' seurat_demoData <- to.Seurat(object = target_demoData, ident = "cell_line", 
#'                              normData = "exprs_norm", forceRaw = FALSE)
#' 
#' @importFrom CreateSeuratObject Seurat
#' @importFrom AddMetaData Seurat
#' 
#' @export 
#' 

to.Seurat <- function(object, ident = NULL, normData = NULL, coordinates = NULL, forceRaw = FALSE){
    
    if (!try(requireNamespace("Seurat", quietly = TRUE))) {
        stop("Please install Seurat from CRAN before converting to a Seurat object")
    }else{
        requireNamespace("Seurat", quietly = TRUE)
    }
    
    
    if(object@featureType == "Probe"){
        stop("Data must be on Target level before converting to a Seurat Object")
    }
    
    if(is.null(normData)){
        stop("normData must not be NULL, please provide name of normalized counts matrix")
    }
    
    if(!normData %in% names(object@assayData)){
        stop(paste0("The normData name \"", normData, "\" is not a valid assay name. Valid names are: ", 
                    paste(names(object@assayData), collapse = ", ")))
    }
    
    normFactor_names <- "normFactors|qFactors|negFactors|hkFactors|hknormFactors"
    
    if(length(grep(pattern = normFactor_names, names(pData(object)))) == 0 & forceRaw == FALSE){
        stop("It is NOT recommended to use Seurat's normalization for GeoMx data. 
             Normalize using GeomxTools::normalize() or set forceRaw to TRUE if you want to continue with Raw data")
    }
    
    
    sequencingMetrics <- c("FileVersion", "SoftwareVersion", "Date", "Plate_ID", "Well", "SeqSetId", "Raw", "Trimmed", 
                           "Stitched", "Aligned", "umiQ30", "rtsQ30", "DeduplicatedReads", "NTC_ID", "NTC")
    
    seuratConvert <- suppressWarnings(Seurat::CreateSeuratObject(counts = object@assayData[[normData]], assay = "GeoMx", 
                                        project = object@experimentData@title))
    seuratConvert <- suppressWarnings(Seurat::AddMetaData(object = seuratConvert, 
                                                  metadata = sData(object)[,!colnames(sData(object)) %in% sequencingMetrics]))
    seuratConvert@assays$GeoMx <- Seurat::AddMetaData(object = seuratConvert@assays$GeoMx, metadata = fData(object))
    
    if(!is.null(ident)){
        if(!ident %in% colnames(seuratConvert@meta.data)){
            stop(paste0("ident \"", ident, "\" not found in GeoMxSet Object"))
        }
        
        Seurat::Idents(seuratConvert) <- seuratConvert[[ident]]
    }
    
    
    seuratConvert@misc <- object@experimentData@other 
    seuratConvert@misc[["sequencingMetrics"]] <- sData(object)[colnames(sData(object)) %in% sequencingMetrics]
    
    if(!is.null(coordinates)){
        xcoord <- coordinates[1]
        ycoord <- coordinates[2]
        
        if(xcoord %in% colnames(seuratConvert@meta.data) & ycoord %in% colnames(seuratConvert@meta.data)){
            coord.df <- data.frame(x=seuratConvert@meta.data[[xcoord]], y=seuratConvert@meta.data[[ycoord]])
            colnames(coord.df) <- coordinates
            seuratConvert@meta.data <- seuratConvert@meta.data[!colnames(seuratConvert@meta.data) %in% coordinates]
        }else{
            if(!xcoord %in% colnames(seuratConvert@meta.data) &
               !ycoord %in% colnames(seuratConvert@meta.data)){
                stop(paste0("xcoord \"", xcoord, "\" and ycoord \"", ycoord, "\" not found in GeoMxSet Object"))
            }
            
            if(!xcoord %in% colnames(seuratConvert@meta.data)){
                stop(paste0("xcoord \"", xcoord, "\" not found in GeoMxSet Object"))
            }
   
            if(!ycoord %in% colnames(seuratConvert@meta.data)){
                stop(paste0("ycoord \"", xcoord, "\" not found in GeoMxSet Object"))
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

#' Convert GeoMxSet Object to SpatialExperiment
#' 
#' @param object GeoMxSet object to convert 
#' @param normData assay containing normalized data
#' @param coordinates X and Y coordinates of each ROI, format: c(X,Y)
#' @param forceRaw should raw data be forced into SpatialExperiment, not recommended
#' 
#' @return SpatialExperiment containing GeoMx data
#' 
#' @examples
#' datadir <- system.file("extdata", "DSP_NGS_Example_Data", package = "GeomxTools")
#' demoData <- readRDS(file.path(datadir, "/demoData.rds"))
#' 
#' target_demoData <- aggregateCounts(demoData)
#' 
#' target_demoData <- normalize(target_demoData, "quant")
#' 
#' seurat_demoData <- to.SpatialExperiment(object = target_demoData, normData = "exprs_norm", 
#'                                         forceRaw = FALSE)
#' 
#' @importFrom SpatialExperiment SpatialExperiment
#' 
#' @export 
#' 

to.SpatialExperiment <- function(object, normData = NULL, coordinates = NULL, forceRaw = FALSE){
    
    if (!try(requireNamespace("SpatialExperiment", quietly = TRUE))) {
        stop("Please install SpatialExperiment from Bioconductor before converting to a SpatialExperiment object")
    }else{
        requireNamespace("SpatialExperiment", quietly = TRUE)
    }
    
    if(object@featureType == "Probe"){
        stop("Data must be on Target level before converting to a SpatialExperiment Object")
    }
    
    if(is.null(normData)){
        stop("normData must not be NULL, please provide name of normalized counts matrix")
    }
    
    if(!normData %in% names(object@assayData)){
        stop(paste0("The normData name \"", normData, "\" is not a valid assay name. Valid names are: ", 
                    paste(names(object@assayData), collapse = ", ")))
    }
    
    normFactor_names <- "normFactors|qFactors|negFactors|hkFactors|hknormFactors"
    
    if(length(grep(pattern = normFactor_names, names(pData(object)))) == 0 & forceRaw == FALSE){
        stop("It is NOT recommended to use Seurat's normalization for GeoMx data. 
             Normalize using GeomxTools::normalize() or set forceRaw to TRUE if you want to continue with Raw data")
    }
    
    sequencingMetrics <- c("FileVersion", "SoftwareVersion", "Date", "Plate_ID", "Well", "SeqSetId", "Raw", "Trimmed", 
                           "Stitched", "Aligned", "umiQ30", "rtsQ30", "DeduplicatedReads", "NTC_ID", "NTC")
    
    if(!is.null(coordinates)){
        xcoord <- coordinates[1]
        ycoord <- coordinates[2]
        
        if(xcoord %in% colnames(sData(object)) & ycoord %in% colnames(sData(object))){
            coord.df <- data.frame(x=sData(object)[[xcoord]], y=sData(object)[[ycoord]])
            colnames(coord.df) <- coordinates
            sequencingMetrics <- c(sequencingMetrics, coordinates)
        }else{
            if(!xcoord %in% colnames(sData(object)) &
               !ycoord %in% colnames(sData(object))){
                stop(paste0("xcoord \"", xcoord, "\" and ycoord \"", ycoord, "\" not found in GeoMxSet Object"))
            }
            
            if(!xcoord %in% colnames(sData(object))){
                stop(paste0("xcoord \"", xcoord, "\" not found in GeoMxSet Object"))
            }
            
            if(!ycoord %in% colnames(sData(object))){
                stop(paste0("ycoord \"", xcoord, "\" not found in GeoMxSet Object"))
            }
        }
        
        rownames(coord.df) <- rownames(sData(object))
        
        coord.df <- as.matrix(coord.df)
        
    }else{
        coord.df <- NULL
    }
    
    spe <- SpatialExperiment::SpatialExperiment(assay = object@assayData[[normData]],
                                                colData = sData(object)[!colnames(sData(object)) %in% 
                                                                            sequencingMetrics],
                                                rowData = fData(object),
                                                spatialCoords = coord.df)
    
    names(spe@assays@data) <- "GeoMx"
    
    spe@metadata <- object@experimentData@other 
    spe@metadata[["sequencingMetrics"]] <- sData(object)[colnames(sData(object)) %in% sequencingMetrics]
    
    return(spe)
}
