#### Specs for GeomxTools:
1. Package shall be downloadable and installable from Bioconductor.
2. Package shall be downloadable and installable from GitHub.
3. The package vignette shall knit into an html file.
4. The knitted vignette shall have no errors in the R code blocks.

#### Specs for readDccFile 
1. Names including target names within the readDccFile returned list elements match the corresponding DCC section labels.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readDccFile.R#L11
2. The number of targets in the returned list matches the number of targets listed in the DCC file.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readDccFile.R#L25
3. The count for each target within the returned list matches the DCC file counts.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readDccFile.R#L33
    
#### Specs for readPkcFile 
1. Column names of PKC files are in correct format.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readPKC.R#L15
2. Names of metadata of PKC files are in correct format.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readPKC.R#L21
3. The number of probes matches the number of probes found in the PKC file.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readPKC.R#L26

#### Specs for readNanoStringGeoMxSet 
1. NanoStringGeoMxSet@assayData$expr dimension labels shall match DCC files input.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readGeoMxSet.R#L35
2. NanoStringGeoMxSet@phenoData@data column names shall match column names in input annotation file and rownames shall match input DCC file names.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readGeoMxSet.R#L44
3. NanoStringGeoMxSet@protocolData@data column names shall match column names in input annotation file and rownames shall match input DCC file names.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readGeoMxSet.R#L65
4. NanoStringGeoMxSet@featureData@data probes shall match probes in the input PKC file.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readGeoMxSet.R#L79
5. Experiment data, testData@experimentData@other, shall contain expected headers.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readGeoMxSet.R#L87
6. The probe counts in testData@assayData$exprs shall match those in the input DCC files.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readGeoMxSet.R#L101
7. Data shall match corresponding data when the same experiment is loaded in DSPDA.     
test:
8. User shall be able to pass optional parameters to read annotation file in with user-specified settings.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readGeoMxSet.R#L136
9. Unidentifiable optional parameters for annotation input reading shall result in an error.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readGeoMxSet.R#L143
10. PKCs are removed if they aren't in the provided config file.
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readGeoMxSet.R#L165
11. Only a single analyte is read in to a GeomMxSet object.
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_readGeoMxSet.R#L209
    
#### Specs for NanoStringGeoMxSet-class
1. The sData row names shall match DCC file names and column names shall match phenoData and protocolData column names.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_accessors.R#L27
2. svarLabels list matches the column names of the sData data frame.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_accessors.R#L35
3. dimLabels shall return a vector with two elements, TargetName and SampleID, by default.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_accessors.R#L43
4. The design slot shall be accessible and NULL by default.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_accessors.R#L52
5. featureType shall be accessible before and after aggregation with default value set to Probe and after aggregation set to Target.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_accessors.R#L59
6. The function countsShiftedByOne shall return a boolean indicating if count shift has been performed on the object.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_accessors.R#L66
7. The function analyte shall return which analyte is read into the object.
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_accessors.R#L76
8. dimLabels replacement method shall allow the replacement of the vector elements within dimLabels.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_replacers.R#L13
9. The design replacement method shall allow the assignment of the design.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_replacers.R#L21
10. The featureType replacement method shall allow the user to replace the feature type of the object.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_replacers.R#L30
    
#### Specs for aggregateCounts
1. All collapsed probe values must be equal to the geoMean of probes for each target.     
test:https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_probeCollapse.R#L59
2. The geomean and geosd of negatives is correct. 
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_probeCollapse.R#L97
3. featureType shall default to Probe and change to Target after the aggregateCounts method is used.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_probeCollapse.R#L117
4. Aggregated GeoMxSet object shall have the same sample dimensions as the probe object and feature dimension shall match the number and names of targets in the feature data target name column. Aggregated object expression matrix dimension labels shall match the aggregated object dimension labels and counts shall by default match the geometric mean of aggregated counts by target from the probe object.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_probeCollapse.R#L124
5. An aggregation function other than geometric mean shall generate an aggregated object with dimension labels matching the default aggregated object and counts aggregated by the specified function for each target.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_probeCollapse.R#L138
6. Other aggregation functions will work. 
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_probeCollapse.R#L155
7. An error shall be displayed if the GeoMxSet object is target-level.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_probeCollapse.R#L174
8. Targets with only one probe shall not have aggregated counts. 
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_probeCollapse.R#L183
9. A warning shall be displayed if no negatives are found in the probe-level object being aggregated.      
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_probeCollapse.R#L189
10. A warning shall be displayed if the GeoMxSet object has no multiprobe targets.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_probeCollapse.R#L194

#### Specs for summarizeNegatives
1. The summarizeNegatives function shall be able to be run multiple times.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_GeoMxSet_probeCollapse.R#L203
2. For each panel, negatives shall be summarized by geometric mean and geometric standard deviation per panel. Summary results shall be saved in sample data.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_shiftCountsOne.R#L33

    
#### Specs for normalize
1. Normalization factors shall be appended to sample data.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_normalization.R#L19
2. Quantile normalization shall allow normalization at any quantile calculated as: AOI quantile / geometric mean of all AOI quantiles.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_normalization.R#L25
3. Quantile normalization shall normalize expression counts as: count/(quantile normalization factor).     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_normalization.R#L51
4. Negative normalization shall calculate normalization factors as: geometric mean of AOI negatives / geometric mean of all AOIs geometric mean of each AOI. One normalization factor shall be calculated per panel for each AOI.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_normalization.R#L66
5. Negative normalized values shall be calculated as: count / negative normalization factor for the corresponding target panel.      
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_normalization.R#L84
6. Housekeeper normalization shall calculate the normalization factor as: geometric mean of housekeepers within AOI / geometric mean of all AOIs geometric mean of housekeepers within AOI.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_normalization.R#L139
7. The user shall be able to do background subtraction by panel. For each panel, negative background shall be calculated as the geometric mean of negatives for each AOI. Background subtracted counts shall be calculated as  counts - corresponding AOI background for the corresponding panel.      
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_normalization.R#L159
8. The user shall be able to do background subtraction with all negatives regardless of panel. Negative background shall be calculated as the geometric mean of negatives for each AOI. Background subtracted counts shall be calculated as  counts - corresponding AOI background.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/0bbd08db0081eb9347e6a859c98d4e363e883ca4/tests/testthat/test_normalization.R#L175
9. An error shall be displayed if no negatives are available to calculate background when performing background subtraction.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/0bbd08db0081eb9347e6a859c98d4e363e883ca4/tests/testthat/test_normalization.R#L208
10. An error shall be displayed if no negatives are available to calculate background when performing negative normalization.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/0bbd08db0081eb9347e6a859c98d4e363e883ca4/tests/testthat/test_normalization.R#L126
11. An error shall be displayed if negative normalization is attempted on a probe-level object.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/0bbd08db0081eb9347e6a859c98d4e363e883ca4/tests/testthat/test_normalization.R#L133
12. An error shall be displayed if background subtraction is attempted on a probe-level object.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/0bbd08db0081eb9347e6a859c98d4e363e883ca4/tests/testthat/test_normalization.R#L161
13. All normalization methods work on protein data.
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_protein.R#L168

    
#### Specs for Quality Control
1. The Technical signal QC shall flag samples that have raw counts below the raw read threshold parameter. The Low Reads flag shall be stored in the protocol data.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L24
2. The Technical signal QC shall flag samples that are below the percent trimmed allowed. The percent trimmed shall be calculated as: 100% * (trimmed reads/raw reads). The LowTrimmed flag shall be stored in the protocol data.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L30
3. The Technical signal QC shall flag samples that are below the percent stitched reads allowed. The percent stitched reads shall be computed as: 100 * (stitched reads / raw reads). The LowStitched flag shall be stored in the protocol data.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L36
4. The Technical signal QC shall flag samples that have percent aligned reads below the percent aligned reads threshold which is the minimum percent of raw reads aligned to target sequences allowed. The percent aligned reads shall be computed as: 100 * (aligned reads / raw reads). The LowAligned flag shall be stored in the protocol data.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L42
5. The Technical signal QC shall flag samples that are below the percent saturation allowed. The percent saturation shall be calculated as: 100% * (1 - deduplicated reads/aligned reads). The LowSaturation flag shall be stored in the protocol data.      
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L48
6. The Negative probe count geomean threshold establishes the level of technical noise expected and flags segments with signal below that level. The LowNegatives QC shall flag samples with geometric mean of negatives below the negative cutoff.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L69
7. The No Template Control (NTC) Count establishes the level at which counts in the NTC will be flagged. The No Template Control (NTC) is used to detect contamination in the library prep. The QC shall flag segments that have higher than the maximum NTC count cutoff.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L82
8. The minimum nuclei count  threshold is the minimum recommended value for the GeoMx run and QC shall flag segments with nuclei count below that level. The LowNuclei flag shall be set for this QC on the object protocol slot.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L95
9. The minimum surface area  threshold is the minimum recommended value for the GeoMx run and QC shall flag segments with area  below that level. The LowArea flag shall be set for this QC on the object protocol slot.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L101
10. The LowNuclei and/or LowArea flag QC shall be bypassed if no area and/or nuclei columns are in the input annotation file.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L109
11. The LowProbeRatio QC shall flag probes with geometric mean expression to target geometric mean expression ratios below the ratio cutoff.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L156
12. The LocalGrubbsOutlier QC shall flag probes that are Grubbs outliers in each AOI when compared to probes from the same target given specified alpha setting.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L163
13. Targets with less than three probes shall have no Grubbs outlier flags.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L196
14. The GlobalGrubbsOutlier QC shall flag probes with percent of LocalGrubbsOutlier flags greater than the percent cutoff across all AOI.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L205
15. LowProbeRatio flags shall match between GeomxTools BioProbeQC and DSPDA BioProbe QC.     
test:
16. GlobalGrubbs flags shall match between GeomxTools BioProbeQC and DSPDA BioProbe QC.     
test:
17. Proportions of local Grubbs' outliers per target shall match proportions per target identified with local Grubbs' outlier testing in DSPDA. Exact probes selected as outliers may vary if more than one probe has the same extreme value. By the definition of the Grubbs' test only one of these values may be excluded and the choice between probes with the same value is arbitrary.      
test:
18. The QC flags from setSegmentQCFlags shall match the QC flags from running setSeqQCFlags, setBackgroundQCFlags, and setGeoMxQCFlags on the same GeoMxSet object with the same settings.      
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_QC.R#L130
19. The QC flags from setSegmentQCFlags shall match the QC flags from running segment Quality Control in DSPDA with the same settings on the same data.     
test:
20. The appropriate segment QC flags shall be added to a protein dataset.
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_protein.R#L103
21. A warning is given if protein data is run through setBioProbeQCFlags.
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_protein.R#L114
    
    
#### Specs for mixedModelDE
1. The function probe p-value shall match the probe p-value from running linear mixed model analysis outside of the function.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_de.R#L65
2. With a non-Windows operating system, parallelizing across multiple cores or clusters shall produce the same results. Multicore parallelization is not available for Windows.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_de.R#L68
3. With a non-Windows operating system, parallelizing across multiple cores or running the function sequentially shall produce the same results. Multicore parallelization is not available for Windows.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_de.R#L74
4. An error shall be produced if model terms are not found in sample data.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_de.R#L100
5. An error shall be produced if grouping variable is not located in the sample data.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_de.R#L115
    
#### Specs for shiftCountsOne
1. The function with useDALogic param set to TRUE shall add one to all zeros in the expression matrix.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_shiftCountsOne.R#L19
2. The function with useDALogic param set to FALSE shall add one count to all counts in the expression matrix.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_shiftCountsOne.R#L27
    
#### Specs for writeNanoStringGeoMxSet
1. The function shall generate a DCC file for each sample from a probe-level GeoMxSet object expression matrix regardless of processing steps previously performed on the object.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_writeNanoStringGeoMxSet.R#L26
2. The function shall not generate DCC files and an error shall be displayed when attempting to write from a target-level GeoMxSet object regardless of processing steps previously performed on the object.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_writeNanoStringGeoMxSet.R#L98
3. GeoMxSet object generated from the created DCC files shall match the GeoMxSet object used to generate the DCC with the writeNanoStringGeoMxSet function with the exception of the NTC sample data column if included.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_writeNanoStringGeoMxSet.R#L46
    
#### Specs for utilty functions
1. The function ngeoMean shall calculate the geometric mean of values. If values <=0 are provided to the function, values shall be shifted by a threshold value to be >0 prior to calculation. Values that are NA shall be ignored in the calculation.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_utils.R#L33
2. The function ngeoSD shall calculate the geometric standard deviation of values. If values <=0 are provided to the function, values shall be shifted by a threshold value to be >0 prior to calculation. Values that are NA shall be ignored in the calculation.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_utils.R#L42
3. The function logtBase shall calculate the log of values by base parameter. If values <=0 are provided to the function, values shall be shifted by a threshold value to be >0 prior to calculation. Values that are NA shall be ignored in the calculation.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_utils.R#L52

#### Specs for GeoMxSet coercions
1. The coercion of a GeoMxSet object shall warn users when coercing non-normalized data but will coerce when forced.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_coercions.R#L18
2. The coercion of a GeoMxSet object shall only occur on target level data.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_coercions.R#L39
3. The coercion of a GeoMxSet object shall only occur when a valid norm data count matrix is provided by the user.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_coercions.R#L48
4. The coercion of a GeoMxSet object shall copy the wanted data from the GeoMxSet object to the correct location in the coerced object.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_coercions.R#L77
5. The coercion of a GeoMxSet object shall warn users when the coordinate column names are not valid.     
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_coercions.R#L211

#### Specs for Protein functions
1. igg.names shall return the expected target names.    
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_protein.R#L130
2. hk.names shall return the expected target names.    
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_protein.R#L139
3. computeNormalizationFactors and normalize calculations match.
test: https://github.com/Nanostring-Biostats/GeomxTools/blob/062b3d94e4be924138ff5315441ad82360eabe67/tests/testthat/test_protein.R#L191
