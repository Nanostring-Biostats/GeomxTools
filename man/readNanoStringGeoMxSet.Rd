\name{readNanoStringGeoMxSet}
\alias{readNanoStringGeoMxSet}
\concept{NanoStringGeoMxSet}

\title{Read 'NanoStringGeoMxSet'}

\description{
  Create an instance of class \code{\linkS4class{NanoStringGeoMxSet}} by reading
  data from NanoString GeoMx Digital Count Conversion (DCC) data.
}

\usage{
readNanoStringGeoMxSet(dccFiles, pkcFiles, phenoDataFile, 
                     phenoDataSheet, phenoDataDccColName = "Sample_ID",
                     phenoDataColPrefix = "", protocolDataColNames = NULL,
                     experimentDataColNames = NULL, ...)
}

\arguments{
  \item{dccFiles}{A character vector containing the paths to the DCC files.}
  \item{pkcFiles}{An optional character string representing the path to the
    corresponding PKC file.}
  \item{phenoDataFile}{An optional character string representing the path
    to the corresponding phenotypic excel data file.}
  \item{phenoDataSheet}{An optional character string representing the excel 
  sheet name containing the phenotypic data.}
  \item{phenoDataDccColName}{Character string identifying unique sample identifier 
  column in \code{phenoDataFile}.}
  \item{phenoDataColPrefix}{An optional prefix to add to the phenoData
    column names to distinguish them from the names of assayData matrices,
    featureData columns, and protocolData columns.}
  \item{protocolDataColNames}{Character list of column names from \code{phenoDataFile} containing data about the experimental protocol or sequencing data.}
  \item{experimentDataColNames}{Character list of column names from \code{phenoDataFile} containing data about the experiment's meta-data.}
  \item{...}{Optional parameters to pass to \code{readxl::read_xlsx} function for annotation read in}
}

\value{
  An instance of the \code{\linkS4class{NanoStringGeoMxSet}} class.
}

\author{Zhi Yang & Nicole Ortogero}

\seealso{
\code{\link{NanoStringGeoMxSet}}
}

\examples{
# Data file paths
datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
dccFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
pkc <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
sampleAnnotationFile <- file.path(datadir, "annotations.xlsx")

dccFileColumn <- "Sample_ID"

dccSet <- readNanoStringGeoMxSet(dccFiles=dccFiles,
                               pkcFiles=pkc,
                               phenoDataFile=sampleAnnotationFile,
                               phenoDataSheet="CW005",
                               phenoDataDccColName=dccFileColumn,
                               protocolDataColNames=c("aoi", "cell_line", 
                                                      "roi_rep", "pool_rep", 
                                                      "slide_rep"),
                               experimentDataColNames="panel", 
                               phenoDataColPrefix="")

# All data
dccSet <- readNanoStringGeoMxSet(dccFiles, pkcFile = pkc,
                                 phenoDataFile = sampleAnnotationFile,
                                 phenoDataSheet="CW005")
varLabels(dccSet)

# All data with phenoData prefix
dccSetPhenoPrefix <- readNanoStringGeoMxSet(dccFiles, 
                                pkcFile = pkc, 
                                phenoDataFile = sampleAnnotationFile,
                                phenoDataSheet="CW005", 
                                phenoDataColPrefix = "PHENO_")
varLabels(dccSetPhenoPrefix)
}

\keyword{file}
\keyword{manip}
