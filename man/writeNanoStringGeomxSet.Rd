\name{writeNanoStringGeomxSet}
\alias{writeNanoStringGeomxSet}
\concept{NanoStringGeomxSet}

\title{write 'NanoStringGeomxSet'}

\description{
  Take an instance of class \code{\linkS4class{NanoStringGeomxSet}} and write
  NanoString GeoMx Digital Count Conversion (DCC) data.
}

\usage{
writeNanoStringGeomxSet(x, dir = getwd())
}

\arguments{
  \item{x}{A NanoStringGeomxSet object.}
  \item{dir}{A directory path to save all the DCC files.}
}


\author{Zhi Yang & Nicole Ortogero}

\seealso{
\code{\link{NanoStringGeomxSet}}
}

\examples{
# Data file paths
datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
dccFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
pkc <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
sampleAnnotationFile <- file.path(datadir, "annotations.xlsx")

dccFileColumn <- "Sample_ID"

dccSet <- readNanoStringGeomxSet(dccFiles=dccFiles,
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
writeNanoStringGeomxSet(dccSet)

}

\keyword{file}
\keyword{manip}
