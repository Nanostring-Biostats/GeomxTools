
# GeoMxTools

## Overview

The GeoMxTools package contains tools for analyzing data from
NanoString GeoMx Digital Spatial Profiler (DSP). It provides functions
to read, quality control (QC) and normalize starting from Nanostring
DCC and PKC files generated from the NanoString GeoMx DSP.

It contains the definition of the NanoStringGeoMxSet class which
inherits from Biobase’s ExpressionSet class and NanoStringRCCSet class.

## Installation

### Download the release version from Bioconductor
<https://bioconductor.org/packages/release/bioc/html/GeomxTools.html>

### Install the release version from Bioconductor
``` r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")

BiocManager::install(version="release")

BiocManager::install("GeomxTools")
```

### Install the development version from GitHub
``` r
install.packages("devtools")
library("devtools")
devtools::install_github("Nanostring-Biostats/GeomxTools", 
                         build_vignettes = TRUE, ref = "dev")
```

## Documentation
To get started, it is recommended to begin with the example workflow:    
[GeoMx Workflow](https://www.bioconductor.org/packages/release/workflows/vignettes/GeoMxWorkflows/inst/doc/GeomxTools_RNA-NGS_Analysis.html)

To learn how to start developing with GeoMxTools, view documentation for the
version of this package installed in your system, start R and enter:

``` r
browseVignettes("GeomxTools")
```

## Branches
The release version on Bioconductor is the stable version.
<https://bioconductor.org/packages/release/bioc/html/GeomxTools.html>

The devel version on Bioconductor is upstream of master on GitHub.
It is under active development and no guarantee is made on usability
at any given time.

The dev branch on GitHub is under active development and no guarantee 
is made on usability at any given time.

## Citation
Ortogero, N.; Yang, Z.; Vitancol, R.; Griswold, M.; Henderson, D. 
GeomxTools: NanoString GeoMx Tools. R Package Version 2.1.1. 
NanoString Technologies Inc.; Seattle, WA 98109, USA. 2021. 

## License
This project is licensed under the [MIT license](LICENSE).