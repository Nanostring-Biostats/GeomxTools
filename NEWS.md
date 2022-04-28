# GeomxTools 3.0.1
* Documentation correction on release (equivalent to 3.1.1 on devel)

# GeomxTools 3.0.0
* No changes from 2.99.2
* Bioconductor release 3.15 version

# GeomxTools 2.99.2
* documentation updates

# GeomxTools 2.99.1
* documentation updates

# GeomxTools 2.1.9
* documentation updates

# GeomxTools 2.1.8
* subset objects in tests & examples to decrease time

# GeomxTools 2.1.7
* coercions to Seurat and SpatialExperiment are S3 coercions 

# GeomxTools 2.1.6
* added updateObject capability 

# GeomxTools 2.1.5

## New features:
* Handle multiple PKC file versions for a single module
    * Only common probes in all versions for each module will be used
    * Default behavior most recent PKC version for resolving probe assignments
    * User can override default with a specified version file name
* Added SystematicName and GeneID from PKC to feature metadata


# GeomxTools 2.1.4
* Add code for grubbs test from deprecated outliers package

# GeomxTools 2.1.3

## New features:
* Allow protein NGS experiment data reading
    * New slot, analyte, added to refer to analyte type

# GeomxTools 2.1.2

## Bug fixes:
* Bug fix in outlier testing

# GeomxTools 2.1.1

## New features:
* Add coercion to Seurat and SpatialExperiment

# GeomxTools 2.1.0

* No changes from 1.1.4

# GeomxTools 2.0.0

* No changes from 1.1.4
* Bioconductor release 3.14 version

# GeomxTools 1.99.4

* No changes from 1.1.4

# GeomxTools 1.1.4

## Revisions:
* Update license

# GeomxTools 1.1.3

## Revisions:
* Allow for multipanel background correction in `subtractBackground()`
* Allow user to re-run `summarizeNegatives()`
* Enable optional parameters in `readNanoStringGeoMxSet()` for annotation file reading
* Added and revised test cases throughout

## Bug fixes:
* Fixed `writeNanoStringGeoMxSet` to work with current DCC file format

# GeomxTools 1.1.2

## Revisions:
* Allow users to use more than one DCC version
* Speed improvements in `setBioProbeQC()` and `aggregateCounts()`

## Bug fixes:
* Fix error in `setBioProbeQC()` with single panel objects
* Fix mixed model output to reference correct p-value
* Fix thresholding in utility functions to keep format matrix format inputs

# GeomxTools 1.1.1

## New features:
* Differential expression with linear mixed model method `mixedModelDE()`

## Revisions:
* Handle multi-panel normalization

## Bug fixes:
* Fix skipping of vectors that don't meet Grubbs requirements
* Fix build warning from knitr update

# GeomxTools 0.99.4 - concomittant development branch version
v0.99.4 includes changes beyond 1.0.0

## New features:
* New slot FeatureType to indicate if data is probe- or target-level
* Segment QC `setSegmentQCFlags()` and probe QC `setBioProbeQC()`
* Count aggregation method `aggregateCounts()`
* Common GeoMx normalizations `normalize()`
* Log and count thresholding methods added to utils

## Revisions:
* Updated `readDccFile()` to expand dcc file versions accepted
* Allow user to  without auto-aggregating counts to target-level
* Probe annotations attached to featureData with readPKCFile

# GeomxTools 1.1.0

* No changes from 1.0.0

# GeomxTools 1.0.0

* Initial release, includes load with automatic aggregation to target-level

## User notes:
* This version was included in Bioconductor release 3.13

## Citation:
Ortogero, N.; Yang, Z.; Vitancol, R.; Griswold, M.; Henderson, D. 
GeomxTools: NanoString GeoMx Tools. R Package Version 1.0.0. 
NanoString Technologies Inc.; Seattle, WA 98109, USA. 2021. 

# GeomxTools 0.99.0

* Package template creation
