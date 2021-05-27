.dccMetadata <-
  list(schema =
         list("Header" =
                data.frame(labelDescription =
                             c("The version of the file",
                               "The version of the software used to create the file",
                               "The date of the sample"),
                           minVersion = numeric_version(c("0.01", "0.01", "0.01")),
                           row.names =
                             c("FileVersion", "SoftwareVersion", "Date"),
                           stringsAsFactors = FALSE),
              "Scan_Attributes" =
                data.frame(labelDescription =
                             c("The sample ID",
                               "The plate ID",
                               "The well ID"),
                           row.names = 
                             c("ID", "Plate_ID", "Well"),
                           minVersion = numeric_version(c(rep("0.01", 3L))),
                           stringsAsFactors = FALSE),
              "NGS_Processing_Attributes" =
                data.frame(labelDescription =
                             c(NA_character_,
                               NA_character_,
                               NA_character_,
                               NA_character_,
                               NA_character_,
                               NA_character_,
                               NA_character_,
                               NA_character_),
                           minVersion = numeric_version(c(rep("0.01", 8L))),
                           row.names =
                             c("SeqSetId", "Raw", "Trimmed", 
                               "Stitched", "Aligned", "umiQ30", "rtsQ30"),
                           stringsAsFactors = FALSE),
              "Code_Summary" =
                data.frame(labelDescription =
                             c(NA_character_, NA_character_),
                           minVersion = numeric_version(c(rep("0.01", 2L))),
                           row.names = c("RTS_ID", "Count"),
                           stringsAsFactors = FALSE)
         )
  )


.dccMetadata[["protocolData"]] <-
  do.call(rbind,
          unname(head(.dccMetadata[["schema"]], 3L)))[, "labelDescription",
                                                      drop = FALSE]

rownames(.dccMetadata[["protocolData"]])[rownames(.dccMetadata[["protocolData"]]) == "ID"] <- "SampleID"


.codeClassMetadata <-
  c("CodeClass,IsControl,Analyte",
    "Endogenous,FALSE,gx|cnv|fusion",
    "Housekeeping,TRUE,gx|fusion",
    "Positive,TRUE,general",
    "Negative,TRUE,general",
    "Binding,TRUE,general",
    "Purification,TRUE,general",
    "Reserved,TRUE,general",
    "SNV_INPUT_CTL,TRUE,SNV",
    "SNV_NEG,TRUE,SNV",
    "SNV_POS,TRUE,SNV",
    "SNV_UDG_CTL,TRUE,SNV",
    "SNV_PCR_CTL,TRUE,SNV",
    "SNV_REF,FALSE,SNV",
    "SNV_VAR,FALSE,SNV",
    "PROTEIN,FALSE,protein",
    "PROTEIN_NEG,TRUE,protein",
    "PROTEIN_CELL_NORM,TRUE,protein",
    "Restriction Site,TRUE,CNV",
    "Invariant,TRUE,CNV")
.codeClassMetadata <-
  utils::read.csv(textConnection(paste0(.codeClassMetadata, collapse = "\n")),
           colClasses = c("character", "logical", "character"),
           stringsAsFactors = FALSE)


.validDccSchema <-
function(x, fileVersion,
         section = c("Header", "Scan_Attributes", "NGS_Processing_Attributes", "Code_Summary"))
{
  section <- match.arg(section)
  schema <- .dccMetadata[["schema"]][[section]]
  expectedNames <- row.names(schema)[schema[,"minVersion"] <= fileVersion]
  if (identical(colnames(x), expectedNames))
    TRUE
  else
    sprintf("<%s> section must contain %s", section,
            paste0("\"", expectedNames, "\"", collapse = ", "))
}


.allNA <- function(x) {
  all(is.na(x))
}

.allTRUE <- function(x) {
  is.logical(x) && !anyNA(x) && all(x)
}

.allFALSE <- function(x) {
  is.logical(x) && !anyNA(x) && !any(x)
}

.allZero <- function(x) {
  is.numeric(x) && !anyNA(x) && identical(range(x), c(0, 0))
}

.validNonNegativeInteger <- function(x) {
  is.integer(x) && !anyNA(x) && min(x) >= 0L
}

.validNonNegativeNumber <- function(x) {
  is.numeric(x) && !anyNA(x) && min(x) >= 0
}

.validPositiveNumber <- function(x) {
  is.numeric(x) && !anyNA(x) && min(x) > 0
}
