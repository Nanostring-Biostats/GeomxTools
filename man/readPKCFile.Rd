\name{readPKCFile}
\alias{readPKCFile}

\title{Read PKC File}

\description{
  Read a NanoString Probe Kit Configuration (PKC) file.
}

\usage{
readPKCFile(file, default_versions=NULL)
}

\arguments{
  \item{file}{A character string containing the path to the PKC file.}
  \item{default_versions}{Optional list of pkc file names to use if more than one version of module pkcs provided.}
}

\value{
  An instance of the \code{\linkS4class{DataFrame}} class containing columns:
  \item{"RTS_ID"}{unique probe ID}
  \item{"TargetName"}{target or gene name}
  \item{"Module"}{PKC name}
  \item{"Negative"}{negative probe}
  \item{\ldots}{additional columns}
}

\author{Zhi Yang & Nicole Ortogero}

\seealso{\code{\link{readNanoStringGeoMxSet}}}

\examples{
datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="GeomxTools")
pkc <- unzip(zipfile = file.path(datadir,  "/pkcs.zip"))
PKCData <- readPKCFile(pkc)
}

\keyword{file}
\keyword{manip}
