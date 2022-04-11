#### Reqs for GeomxTools:
- The GeomxTools package shall be downloaded and installed from Bioconductor or from GitHub without GeomxTools-specific installation errors.
- The package vignette shall knit without error.

Specifications:

#### Reqs for readDccFile 
- The function shall read in a Digital Count Conversion (DCC) file for one AOI 
- The function shall return a list containing: 
  - Header - version numbers and date
  - Scan Attributes - plate and well information
  - NGS Processing Attributes - sequencing information
  - Code Summary - counts per RTS_ID 
      
Specifications:
    
#### Reqs for readPkcFile 
- The function shall read in one or more Probe Kit Configuration (PKC) files 
- The function shall return a dataframe with 6 columns and Metadata:
  - Column Names:
    - RTS_ID
    - Target (Gene Name)
    - Module (PKC file name)
    - CodeClass (Gene Type: Negative, Control, Endogenous, etc)
    - ProbeID
    - Negative (Negative Probe Boolean)
  - Metadata:
    - PKC File Name(s)
    - PKC File Version(s)
    - PKC File Date(s)
    - Analyte Type
    - Minimum AOI Area Recommendation(s)
    - Minimum Nuclei Count Recommendation(s)
      
Specifications:
    
#### Reqs for readNanoStringGeoMxSet 
- The function shall return a NanoStringGeoMxSet object with the following attributes:
  - Assay Data - gene x AOI count matrix
  - Pheno Data - AOI segment properties 
  - Protocol Data - AOI annotations
  - Feature Data - read PKC output (probe annotations)
  - Experiment Data - lab, PI, publication data, PKC metadata, and any user-defined columns
  
Specifications:
    
#### Reqs for NanoStringGeoMxSet-class
- There shall be functions to access parts of NanoStringGeoMxSet by users
  - sData
  - svarLabels
  - dimLabels
  - design
  - featureType
  - countsShiftedByOne
- There shall be functions to replace parts of NanoStringGeoMxSet by users
  - dimLabels
  - design
  - featureType
  
Specifications:
    
#### Reqs for aggregateCounts
- The function shall collapse multi-probe targets using geometric mean by default or a user-specified aggregation function.

Specifications:
    
#### Reqs for summarizeNegatives
- The function shall add negative geometric mean, negative geometric standard deviation, and any additional negative aggregations specified by the user to sample annotations.

Specifications:
    
#### Reqs for normalize
- The method shall add normalized data as an assayDataElement matrix with quantile normalization, negative normalization, housekeeper normalization, or background subtraction.

Specifications:
    
#### Reqs for Quality Control
- GeomxTools shall provide the ability to perform Segment QC on the data by performing Technical Signal quality Control which assesses the quality sequencing of each segment.  
- GeomxTools shall provide the ability to perform Segment QC on the data by performing Technical background quality Control which flags the samples based on the negative controls.  
- GeomxTools shall provide the ability to perform Segment QC on the data based on minimum nuclei and surface area count cutoffs. 
  - These QCs may be safely bypassed if nuclei or area data is not available.
- The setSegmentQCFlags method shall be equivalent to running all segment QC methods (setSeqQCFlags, setBackgroundQCFlags, setGeoMxQCFlags); QC results shall match segement Qualty Control in DSPDA performed with the same settings.
- Biological probe QC shall be available in GeomxTools. This shall allow probes that appear to be outliers in the data to be flagged and results shall be similar to BioProbe Qualty Control in DSPDA performed with the same settings.

Specifications:
    
#### Reqs for mixedModelDE
- The function shall perform linear mixed model different expression analysis on expression data.

Specifications:
    
#### Reqs for shiftCountsOne
- The function shall impute counts by adding one.

Specifications:
    
#### Reqs for writeNanoStringGeoMxSet
- The function shall save expression data from a GeoMxSet object as DCC file(s).

Specifications:
    
#### Reqs for utilty functions
- The function ngeoMean shall perform geometric mean transformations.
- The function ngeoSD shall perform geometric standard deviation transformations.
- The function logtBase shall perform log transformations.

Specifications:

#### Reqs for coercions
- The function to.Seurat shall copy neccesary data to a Seurat object.
- The function to.SpatialExperiment shall copy neccesary data to a SpatialExperiment object.

Specifications:
