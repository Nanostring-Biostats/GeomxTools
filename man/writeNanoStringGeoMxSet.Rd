\name{writeNanoStringGeoMxSet}
\alias{writeNanoStringGeoMxSet}
\concept{NanoStringGeoMxSet}

\title{write 'NanoStringGeoMxSet'}

\description{
  Take an instance of class \code{\linkS4class{NanoStringGeoMxSet}} and write
  NanoString GeoMx Digital Count Conversion (DCC) data.
}

\usage{
writeNanoStringGeoMxSet(x, dir = getwd())
}

\arguments{
  \item{x}{A NanoStringGeoMxSet object.}
  \item{dir}{A directory path to save all the DCC files.}
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
writeNanoStringGeoMxSet(dccSet)

}

\keyword{file}
\keyword{manip}
