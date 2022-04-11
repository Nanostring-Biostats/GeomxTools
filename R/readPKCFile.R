readPKCFile <-
function(file, default_pkc_vers=NULL)
{
  pkc_json_list <- lapply(file, function(pkc_file) {rjson::fromJSON(file = pkc_file)})
  pkc_names <- extract_pkc_names(file)
  names(pkc_json_list) <- pkc_names
  pkc_modules <- basename(unlist(lapply(pkc_names, sub, pattern="_[^_]+$", replacement="")))
  names(pkc_modules) <- pkc_names
  # Extract header
  header <- list(PKCFileName = sapply(pkc_json_list, function(list) list[["Name"]]),
                 PKCModule = pkc_modules,
                 PKCFileVersion = sapply(pkc_json_list, function(list) list[["Version"]]),
                 PKCFileDate = sapply(pkc_json_list, function(list) list[["Date"]]),
                 AnalyteType = sapply(pkc_json_list, function(list) list[["AnalyteType"]]),
                 MinArea = sapply(pkc_json_list, function(list) list[["MinArea"]]),
                 MinNuclei = sapply(pkc_json_list, function(list) list[["MinNuclei"]])  
                 )

  # Check for multiple versions of pkc
  multi_idx <- duplicated(header[["PKCModule"]])
  multi_mods <- unique(header[["PKCModule"]][multi_idx])
  if (length(multi_mods) < 1) {
    if (!is.null(default_pkc_vers)) {
      warning("Only one version found per PKC module. ",
              "No PKCs need to be combined. ",
              "Therefore, no default PKC versions will be used.")
    }
  } else {
    use_pkc_names <- lapply(multi_mods, function(mod) {
        mod_idx <- header[["PKCModule"]] == mod
        max_vers <- as.numeric(as.character(max(as.numeric_version(
          header[["PKCFileVersion"]][mod_idx]))))
        max_name <- names(header[["PKCFileVersion"]][
          header[["PKCFileVersion"]] == max_vers])
        return(max_name)
      })
    names(use_pkc_names) <- multi_mods
    if (!is.null(default_pkc_vers)) {
      default_names <- extract_pkc_names(default_pkc_vers)
      default_mods <- extract_pkc_modules(default_pkc_vers)
      dup_defaults <- default_names[duplicated(default_mods) | 
        duplicated(default_mods, fromLast=TRUE)]
      if (!all(default_names %in% names(header[["PKCFileName"]]))) {
        removed_pkcs <- 
          default_pkc_vers[!default_names %in% names(header[["PKCFileName"]])]
        stop("Could not match all default PKC versions with a PKC file name. ", 
             "Check default file names match exactly to a PKC file name.\n",
             paste0("Unmatched default PKC versions: ", removed_pkcs))
      } else if (length(dup_defaults) > 0) {
        stop("There should only be one default PKC version per module. ", 
             "Ensure only one version per module in default PKCs list.\n",
             "Multiple default PKC version conflicts: ", 
             paste(dup_defaults, collapse=", "))
      } else {
        use_pkc_names[default_mods] <- default_names
      }
    }
  }

  rtsid_lookup_df <- generate_pkc_lookup(pkc_json_list)
  # create negative column 
  rtsid_lookup_df$Negative <- grepl("Negative", rtsid_lookup_df$CodeClass)
  rtsid_lookup_df$RTS_ID <- gsub("RNA", "RTS00", rtsid_lookup_df[["RTS_ID"]])
  # Coerce output to DataFrame
  rtsid_lookup_df <- S4Vectors::DataFrame(rtsid_lookup_df)

  if (length(multi_mods) > 0) {
    for (mod in names(use_pkc_names)) {
      mod_vers <- names(header[["PKCModule"]])[header[["PKCModule"]] == mod]
      mod_lookup <- rtsid_lookup_df[rtsid_lookup_df$Module %in% mod_vers, ]
      mod_tab <- table(mod_lookup$RTS_ID)
      remove_rts <- names(mod_tab[mod_tab != length(mod_vers)])
      if (length(remove_rts) > 0) {
                warning("The following probes were removed from analysis",
                " as they were not found in all PKC module versions used.\n",
                paste(capture.output(print(
                  subset(rtsid_lookup_df, subset=RTS_ID %in% remove_rts))),
                  collapse = "\n"))
        rtsid_lookup_df <- 
          subset(rtsid_lookup_df, subset=!RTS_ID %in% remove_rts)
      }
      remove_vers <- mod_vers[mod_vers != use_pkc_names[mod]]
      rtsid_lookup_df <- 
        subset(rtsid_lookup_df, subset=!Module %in% remove_vers)
      warning("The following PKC versions were removed from analysis",
        " as they were overridden by a newer PKC version or",
        " were overridden by a user-defined default PKC version.\n",
        paste(remove_vers, collapse = ", "))
      header <- lapply(header, function(elem) {
        elem[!names(elem) %in% remove_vers]})
    }
  }

  S4Vectors::metadata(rtsid_lookup_df) <- header

  return(rtsid_lookup_df)
}


generate_pkc_lookup <- function(jsons_vec) {
  lookup_df <- data.frame(RTS_ID=character(), 
                          Target=character(), 
                          Module=character(), 
                          CodeClass=character(), 
                          ProbeID=character(), 
                          GeneID=character(), 
                          SystematicName=character(), 
                          stringsAsFactors=FALSE)
  for (curr_idx in seq_len(length(jsons_vec))) {
    curr_module <- names(jsons_vec)[curr_idx]
    curr_json <- jsons_vec[[curr_idx]]
    for (targ in curr_json[["Targets"]]) {
      curr_targ <- targ[["DisplayName"]]
      curr_code_class <- gsub("\\d+$", "", targ[["CodeClass"]])
      for (prb in targ[["Probes"]]) {
        if(curr_json[["AnalyteType"]] == "Protein"){
          curr_RTS_ID <- targ$RTS_ID
        }else{
          curr_RTS_ID <- prb$RTS_ID
        }
        curr_probe_ID <- prb$ProbeID
        curr_gene_ID <- 
          paste(prb$GeneID, collapse = ", ")
        if (length(prb$GeneID) < 1) {
          curr_gene_ID <- NA
        }
        curr_syst_name <- paste(prb$SystematicName, collapse = ", ")
        lookup_df[nrow(lookup_df) + 1, ] <- 
          list(curr_RTS_ID, curr_targ, curr_module, curr_code_class, 
               curr_probe_ID, curr_gene_ID, curr_syst_name)
      }
    }
  }
  return(lookup_df)
}

generate_pkc_targ_notes <- function(jsons_vec, lookup_tab) {
  # Create non-duplicated map from target to pool and codeclass
  sub_lookup <- unique(lookup_tab[, names(lookup_tab) != "RTS_ID"])
  #rownames(sub_lookup) <- sub_lookup[["Target"]]
  notes_df <- 
    data.frame(TargetName=sub_lookup[["Target"]],
               HUGOSymbol=sub_lookup[["Target"]],
               TargetGroup=rep("All Probes", length(rownames(sub_lookup))),
               AnalyteType=rep("RNA", nrow(sub_lookup)),
               CodeClass=sub_lookup[, "CodeClass"],
               Pooling=sub_lookup[, "Module"],
               stringsAsFactors=FALSE)
  for (curr_idx in seq_len(length(jsons_vec))) {
    curr_module <- names(jsons_vec)[curr_idx]
    curr_json <- jsons_vec[[curr_idx]]
    if(length(curr_json[["ProbeGroups"]]) > 0) {
      for (prb_group in curr_json[["ProbeGroups"]]) {
        curr_group <- prb_group[["Name"]]
        for (targ in prb_group[["Targets"]]) {
          notes_df[notes_df$TargetName == targ, "TargetGroup"] <-
            paste(notes_df[notes_df$TargetName == targ, "TargetGroup"], 
                  curr_group, sep=";")
        }
      }
    }
  }
  
  return(notes_df)
}

extract_pkc_names <- function(pkc_files) {
  pkc_names <- 
    unlist(lapply(pkc_files, function(pkc_file) {
      base_pkc_name <- gsub(".pkc", "", trimws(basename(pkc_file)))
      return(base_pkc_name)
    }))
  return(pkc_names)
}

extract_pkc_versions <- function(pkc_files) {
  pkc_files <- extract_pkc_names(pkc_files)
  pkc_vers <- unlist(lapply(pkc_files, sub, pattern="^.*_v", replacement=""))
  return(pkc_vers)
}

extract_pkc_modules <- function(pkc_files) {
  pkc_files <- extract_pkc_names(pkc_files)
  pkc_vers <- unlist(lapply(pkc_files, sub, pattern="_[^_]+$", replacement=""))
  return(pkc_vers)
}
