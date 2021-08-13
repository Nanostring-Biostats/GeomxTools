\name{readDccFile}
\alias{readDccFile}

\title{Read DCC File}

\description{
  Read a NanoString GeoMx Digital Count Conversion (DCC) file.
}

\usage{
readDccFile(file)
}

\arguments{
  \item{file}{A character string containing the path to the DCC file.}
}

\value{
  A list object with two elements:
  \item{\code{"Header"}}{a \code{data.frame} object containing the protocol and sequencing
    information.}
  \item{\code{"Code_Summary"}}{a \code{data.frame} object containing the
    target probe counts.}
}

\author{Zhi Yang & Nicole Ortogero}

\seealso{\code{\link{readNanoStringGeoMxSet}}}

\examples{
datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
dccFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
dccData <- sapply(dccFiles, readDccFile, simplify = FALSE)
}

\keyword{file}
\keyword{manip}
