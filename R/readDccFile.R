readDccFile <-
function(file)
{
  # Read data from Reporter Code Count (RCC) file
  lines <- trimws(readLines(file))

  # Split data by tags
  tags <- names(.dccMetadata[["schema"]])
  output <- sapply(tags, function(tag)
  {
    bounds <- charmatch(sprintf(c("<%s>", "</%s>"), tag), lines)
    if (anyNA(bounds) || bounds[1L] + 1L >= bounds[2L])
      lines[integer(0)]
    else
      lines[(bounds[1L] + 1L):(bounds[2L] - 1L)]
  }, simplify = FALSE)

  # Convert single row attributes to data.table objects
  for (tag in c("Header", "Scan_Attributes", "NGS_Processing_Attributes")) {
    while (length(bad <- grep(",", output[[tag]], invert = TRUE)) > 0L) {
      bad <- bad[1L]
      if (bad == 1L)
        stop(sprintf("%s section has malformed first line", tag))
      fixed <- output[[tag]]
      fixed[bad - 1L] <- sprintf("%s %s", fixed[bad - 1L], fixed[bad])
      output[[tag]] <- fixed[- bad]
    }
    output[[tag]] <- strsplit(output[[tag]], split = ",")
    output[[tag]] <-
      structure(lapply(output[[tag]],
                       function(x) if (length(x) == 1L) "" else x[2L]),
                names = lapply(output[[tag]], `[`, 1L),
                class = "data.frame",
                row.names = basename(file))
  }

  # Coerce numeric_version information in Header
  metadata <- .dccMetadata[["schema"]][["Header"]]
  cols <- c("FileVersion", "SoftwareVersion")
  if (!(all(cols %in% colnames(output[["Header"]]))))
    stop("Header section must contain \"FileVersion\" and \"SoftwareVersion\"")
  output[["Header"]][, cols] <- lapply(output[["Header"]][, cols],
                                       numeric_version)

  # Extract FileVersion for internal checks
  fileVersion <- output[["Header"]][1L, "FileVersion"]
  if (!(numeric_version(fileVersion) %in% numeric_version(c("0.01"))))
    stop("\"FileVersion\" in Header section must be 0.01")

  # Check single row attributes
  for (section in c("Header", "Scan_Attributes", "NGS_Processing_Attributes")) {
    valid <- .validDccSchema(output[[section]], fileVersion, section)
    if (!isTRUE(valid))
      stop(valid)
  }

  # Coerce date data type in Sample_Attributes
  output[["Header"]][["Date"]] <-
    as.Date(output[["Header"]][["Date"]], format = "%Y-%m-%d")

  # Coerce numeric data in Lane Attributes
  cols <- c( "Raw", "Trimmed", "Stitched", "Aligned", "umiQ30", "rtsQ30" )
  output[["NGS_Processing_Attributes"]][, cols] <-
    lapply(output[["NGS_Processing_Attributes"]][, cols], as.integer)

  # Coerce the column name ID to be SampleID
  names(output[["Scan_Attributes"]])[names(output[["Scan_Attributes"]]) == "ID"] <- "SampleID"
  
  # Convert Code_Summary to data.frame object
  output[["Code_Summary"]] <- paste0("RTS_ID,Count\n", 
                                     paste(output[["Code_Summary"]], collapse = "\n"))
  output[["Code_Summary"]] <-
    utils::read.csv(textConnection(output[["Code_Summary"]]),
             colClasses = c(RTS_ID = "character", Count = "numeric"))
  output[["Code_Summary"]][["Count"]] <-
    as.integer(round(output[["Code_Summary"]][["Count"]]))
  rn <- output[["Code_Summary"]][["RTS_ID"]]
  if ((ndups <- anyDuplicated(rn)) > 0L) {
    warning(sprintf("removed %d rows from \"Code_Summary\" due to duplicate rownames",
                    ndups))
    ok <- which(!duplicated(rn, fromLast = FALSE) &
                !duplicated(rn, fromLast = TRUE))
    rn <- rn[ok]
    output[["Code_Summary"]] <- output[["Code_Summary"]][ok, , drop = FALSE]
  }
  rownames(output[["Code_Summary"]]) <- rn

  return( output )
}
