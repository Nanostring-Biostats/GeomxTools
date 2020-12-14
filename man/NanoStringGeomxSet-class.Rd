\name{NanoStringGeomxSet-class}
\docType{class}

\alias{class:NanoStringGeomxSet}
\alias{NanoStringGeomxSet-class}

% constructors
\alias{NanoStringGeomxSet}
\alias{NanoStringGeomxSet,environment-method}
\alias{NanoStringGeomxSet,ExpressionSet-method}
\alias{NanoStringGeomxSet,matrix-method}
\alias{NanoStringGeomxSet,missing-method}
\alias{NanoStringGeomxSet,NanoStringGeomxSet-method}

% coercion
\alias{coerce,ExpressionSet,NanoStringGeomxSet-method}

% show
\alias{show,NanoStringGeomxSet-method}

% accessing
\alias{sData}
\alias{sData,NanoStringGeomxSet-method}
\alias{svarLabels}
\alias{svarLabels,NanoStringGeomxSet-method}
\alias{dimLabels,NanoStringGeomxSet-method}
\alias{dimLabels<-,NanoStringGeomxSet,character-method}
\alias{signatures}
\alias{signatures,NanoStringGeomxSet-method}
\alias{signatures<-}
\alias{signatures<-,NanoStringGeomxSet,SignatureSet-method}
\alias{signatureScores}
\alias{signatureScores,NanoStringGeomxSet-method}
\alias{design,NanoStringGeomxSet-method}
\alias{design<-,NanoStringGeomxSet,ANY-method}
\alias{design<-,NanoStringGeomxSet,NULL-method}
\alias{design<-,NanoStringGeomxSet,formula-method}

\title{Class to Contain NanoString Spatial Expression Level Assays}

\description{
  The \code{NanoStringGeomxSet} class extends the
  \code{\linkS4class{ExpressionSet}} class for NanoString GeoMx Digital Count Conversion
  (DCC) data.
}

\usage{
NanoStringGeomxSet(assayData,
                 phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                 featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                 experimentData=MIAME(),
                 annotation=character(),
                 protocolData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                 dimLabels=c("GeneName", "SampleID"),
                 signatures=SignatureSet(),
                 design=NULL,
                 \ldots)
}

\arguments{
  \item{assayData}{A \code{matrix} or \code{environment} containing the DCCs.}
  \item{phenoData}{An \code{\link{AnnotatedDataFrame}} containing the phenotypic
    data of areas of interest.}
  \item{featureData}{An \code{\link{AnnotatedDataFrame}} containing gene information; gene name, accession number, functional groups, etc.}
  \item{experimentData}{An optional \code{\link{MIAME}} instance with meta-data
    about the experiment.}
  \item{annotation}{A \code{character} string for the PKC file(s).}
  \item{protocolData}{An \code{\link{AnnotatedDataFrame}} containing meta-data 
  about the protocol and sequencing; columns could include \code{"FileVersion"},
  \code{"SoftwareVersion"}, \code{"Date"}, \code{"Plate_ID"}, \code{"Well"}, 
  \code{"SeqSetId"}, \code{"trimGaloreOpts"}, \code{"flash2Opts"}, 
  \code{"umiExtractOpts"}, \code{"boxtie2Opts"}, \code{"Raw"}, \code{"Trimmed"},
  \code{"Stitched"}, \code{"Aligned"}, \code{"umiQ30"}, \code{"rtsQ30"}.}
  \item{dimLabels}{A \code{character} vector of length 2 that provides the
    column names to use as labels for the features and samples respectively
    in the \code{autoplot} method.}
  \item{signatures}{An optional \code{\link{SignatureSet}} object containing
    signature definitions.}
  \item{design}{An optional one-sided formula representing the experimental
    design based on columns from \code{\link[Biobase]{phenoData}}}
  \item{\ldots}{Additional arguments for \code{\link{ExpressionSet}}.}
}

\section{Accessing}{
  In addition to the standard \code{\linkS4class{ExpressionSet}} accessor
  methods, \code{NanoStringGeomxSet} objects have the following:
  \describe{
    \item{}{
      \code{sData(object)}: extracts the \code{data.frame} containing the
      sample data, \code{cbind(pData(object), pData(protocolData(object)))}.
    }
    \item{}{
      \code{svarLabels(object)}: extracts the sample data column names,
      \code{c(varLabels(object), varLabels(protocolData(object)))}.
    }
    \item{}{
      \code{dimLabels(object)}: extracts the column names to use as labels
      for the features and samples.
    }
    \item{}{
      \code{dimLabels(object) <- value}: replaces the \code{dimLabels} of
      the \code{object}.
    }
    \item{}{
      \code{signatures(object)}: extracts the \code{\link{SignatureSet}}
      of the \code{object}.
    }
    \item{}{
      \code{signatures(object) <- value}: replaces the
      \code{\link{SignatureSet}} of the \code{object}.
    }
    \item{}{
      \code{signatureScores(object, elt="exprs")}: extracts the matrix
      of computed signature scores.
    }
    \item{}{
      \code{design(object)}: extracts the one-sided formula representing
      the experimental design based on columns from
      \code{\link[Biobase]{phenoData}}.
    }
    \item{}{
      \code{design(object) <- value}: replaces the one-sided formula
      representing the experimental design based on columns from
      \code{\link[Biobase]{phenoData}}.
    }
  }
}

\author{Zhi Yang & Nicole Ortegero}

\seealso{
  \code{\link{readNanoStringGeomxSet}},
  \code{\linkS4class{ExpressionSet}}
}

\examples{
# Create NanoStringGeomxSet from data files
datadir <- system.file("extdata", "DSP_NGS_Example_Data",
                       package="NanoStringNCTools")
dccFiles <- dir(datadir, pattern=".dcc$", full.names=TRUE)
pkc <- dir(datadir, pattern=".pkc$", full.names=TRUE)
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

# Accessing sample data and column names
head(sData(dccSet))
svarLabels(dccSet)

# Accessing number of samples and features
dim(dccSet)
}

\keyword{classes}
\keyword{methods}