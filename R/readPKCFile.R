readPKCFile <-
function(file)
{
  pkc_json_list <- lapply(file, function(pkc_file) {rjson::fromJSON(file = pkc_file)})
  pkc_names <- 
    unlist(lapply( file, 
                   function(file) {
                     base_pkc_name <- gsub(".pkc", "", trimws(basename(file)))
                     return(base_pkc_name)
                   }))
  names(pkc_json_list) <- pkc_names
  rtsid_lookup_df <- generate_pkc_lookup(pkc_json_list)
  # create negative column 
  rtsid_lookup_df$Negative <- grepl("Negative", rtsid_lookup_df$CodeClass)
  rtsid_lookup_df$RTS_ID <- gsub("RNA", "RTS00", rtsid_lookup_df[["RTS_ID"]])
  # Coerce output to DataFrame
  rtsid_lookup_df <- S4Vectors::DataFrame(rtsid_lookup_df)
  
  # Extract header
  header <- list(PKCFileName = sapply(pkc_json_list, function(list) list[["Name"]]),
                 PKCFileVersion = sapply(pkc_json_list, function(list) list[["Version"]]),
                 PKCFileDate = sapply(pkc_json_list, function(list) list[["Date"]]),
                 AnalyteType = sapply(pkc_json_list, function(list) list[["AnalyteType"]]),
                 MinArea = sapply(pkc_json_list, function(list) list[["MinArea"]]),
                 MinNuclei = sapply(pkc_json_list, function(list) list[["MinNuclei"]])  
                 )
  
  S4Vectors::metadata(rtsid_lookup_df) <- header
  
  return(rtsid_lookup_df)
}


generate_pkc_lookup <- function(jsons_vec) {
  lookup_df <- data.frame(RTS_ID=character(), 
                          Target=character(), 
                          Module=character(), 
                          CodeClass=character(), 
                          ProbeID=character(),
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
        lookup_df[nrow(lookup_df) + 1, ] <- 
          list(curr_RTS_ID, curr_targ, curr_module, curr_code_class, curr_probe_ID)
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

