writeNanoStringGeoMxSet <- function(x, dir = getwd()) {
  stopifnot(is(x, "NanoStringGeoMxSet"))
  validObject(x)
  if (!dir.exists(dir)) 
    dir.create(dir)
  features <- pData(featureData(x))[, c("RTS_ID")]
  header <- "<Header>\nFileVersion,%s\nSoftwareVersion,%s\nDate,%s\n</Header>\n"
  scanAttr <- paste0("<Scan_Attributes>\nID,%s\nPlate_ID,%s\nWell,%s\n</Scan_Attributes>\n")
  laneAttr <- paste0("<NGS_Processing_Attributes>\nSeqSetId,%s\ntamperedIni,%s\ntrimGaloreOpts,%s\n", 
                     "flash2Opts,%s\numiExtractOpts,%s\nbowtie2Opts,%s\n", 
                     "umiDedupOpts,%s\nRaw,%d\nTrimmed,%d\n", 
                     "Stitched,%d\nAligned,%d\numiQ30,%.4f\n", 
                     "rtsQ30,%.4f\n</NGS_Processing_Attributes>\n")
  for (i in seq_len(dim(x)[["Samples"]])) {
    protocolRow <- pData(protocolData(x))[i, ]
    fname <- file.path(dir, sampleNames(x)[i])
    if (file.exists(fname)) {
      file.remove(fname)
    }
    con <- file(fname, open = "a")
    writeLines(sprintf(header, protocolRow[["FileVersion"]], protocolRow[["SoftwareVersion"]],
                       protocolRow[["Date"]]), 
                 con)
    writeLines(sprintf(scanAttr, protocolRow[["SampleID"]], protocolRow[["Plate_ID"]], 
                       protocolRow[["Well"]]), con)

    writeLines(sprintf(laneAttr, protocolRow[["SeqSetId"]], protocolRow[["tamperedIni"]], 
                       protocolRow[["trimGaloreOpts"]], protocolRow[["flash2Opts"]], 
                       protocolRow[["umiExtractOpts"]], protocolRow[["bowtie2Opts"]], 
                       protocolRow[["umiDedupOpts"]], protocolRow[["Raw"]],
                       protocolRow[["Trimmed"]], protocolRow[["Stitched"]], 
                       protocolRow[["Aligned"]], protocolRow[["umiQ30"]],
                       protocolRow[["rtsQ30"]]), 
               con)
    writeLines("<Code_Summary>", con)
    write.table(cbind(features[which(exprs(x)[, i]>0)], 
                      Count = exprs(x)[which(exprs(x)[, i]>0), i]), file = con, quote = FALSE, 
                sep = ",", row.names = FALSE, col.names = FALSE)
    writeLines("</Code_Summary>\n", con)
    writeLines("SOMEHASH100000000000", con)
    close(con)
  }
  invisible(file.path(dir, sampleNames(x)))
}
