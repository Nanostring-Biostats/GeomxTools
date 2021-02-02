
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GeomxTools

## Overview

The GeomxTools is a package that contains tools for analyzing data from
NanoString GeoMx Digital Spatial Profiler (DSP). It provides functions
to read, perform quality control (QC) and normalization on Nanostring
DCC and PKC files generated from the NanoString GeoMx DSP.

It contains the definition of the NanoStringGeomxSet class which
inherits from Biobase’s ExpressionSet class and NanoStringRCCSet class.

## Installation

You can download the the package from bioconductor from this link
<https://bioconductor.org/packages/devel/bioc/html/GeomxTools.html>

``` r
# Install from Bioconductor
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install(version="devel")

BiocManager::install("GeomxTools")

# Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("Nanostring-Biostats/GeomxTools", 
                         build_vignettes = TRUE)
```

## Documentation

To learn how to start using GeomxTools, view documentation for the
version of this package installed in your system, start R and enter:

``` r
browseVignettes("GeomxTools")
```
