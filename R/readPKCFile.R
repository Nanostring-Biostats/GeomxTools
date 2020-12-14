readPKCFile <-
function(file)
{
  pkc_json_list <- lapply(file, function(pkc_file) {fromJSON(file = pkc_file)})
  pkc_names <- 
    unlist(lapply( file, 
                   function(file) {
                     base_pkc_name <- gsub(".pkc", "", trimws(basename(file)))
                     return(base_pkc_name)
                   }))
  names(pkc_json_list) <- pkc_names
  rnaid_lookup_df <- generate_pkc_lookup(pkc_json_list)
  target_notes <- generate_pkc_targ_notes(pkc_json_list, rnaid_lookup_df)
  # create negative column 
  rnaid_lookup_df$Negative <- 
    rnaid_lookup_df$Gene %in% 
    target_notes[grep("Negative", target_notes$Codeclass), "TargetName"]
  # Remove codeclass column, only needed for target notes generation
  rnaid_lookup_df <- 
    rnaid_lookup_df[, !colnames(rnaid_lookup_df) %in% c("Codeclass", "PoolNum")]
  # change "RTS00" to "RNA" (only check the first characters, also need to change from RNA to RTS00)
  rnaid_lookup_df$RTS_ID <- gsub("RTS00", "RNA", rnaid_lookup_df[["RTS_ID"]])
  # Coerce output to DataFrame
  rnaid_lookup_df <- DataFrame(rnaid_lookup_df)
  
  # Extract header
  header <- list(PKCFileName = sapply(pkc_json_list, function(list) list[["Name"]]),
                 PKCFileVersion = sapply(pkc_json_list, function(list) list[["Version"]]),
                 PKCFileDate = sapply(pkc_json_list, function(list) list[["Date"]]),
                 AnalyteType = sapply(pkc_json_list, function(list) list[["AnalyteType"]]),
                 MinArea = sapply(pkc_json_list, function(list) list[["MinArea"]]),
                 MinNuclei = sapply(pkc_json_list, function(list) list[["MinNuclei"]])  
                 )
  
  metadata(rnaid_lookup_df) <- header
  
  return(rnaid_lookup_df)
}


generate_pkc_lookup <- function(jsons_vec) {
  lookup_df <- data.frame(RTS_ID=character(), 
                          Gene=character(), 
                          Module=character(), 
                          Codeclass=character(),
                          stringsAsFactors=FALSE)
  for (curr_idx in 1:length(jsons_vec)) {
    curr_module <- names(jsons_vec)[curr_idx]
    curr_json <- jsons_vec[[curr_idx]]
    for (targ in curr_json[["Targets"]]) {
      curr_gene <- targ[["DisplayName"]]
      curr_code_class <- targ[["CodeClass"]]
      for (prb in targ[["Probes"]]) {
        curr_RTS_ID <- prb$RTS_ID
        lookup_df[nrow(lookup_df) + 1, ] <- 
          list(curr_RTS_ID, curr_gene, curr_module, curr_code_class)
      }
    }
  }
  return(lookup_df)
}

generate_pkc_targ_notes <- function(jsons_vec, lookup_tab) {
  # Create non-duplicated map from gene to pool and codeclass
  sub_lookup <- unique(lookup_tab[, names(lookup_tab) != "RTS_ID"])
  #rownames(sub_lookup) <- sub_lookup[["Gene"]]
  notes_df <- 
    data.frame(TargetName=sub_lookup[["Gene"]],
               HUGOSymbol=sub_lookup[["Gene"]],
               TargetGroup=rep("All Probes", length(rownames(sub_lookup))),
               AnalyteType=rep("RNA", nrow(sub_lookup)),
               Codeclass=sub_lookup[, "Codeclass"],
               Pooling=sub_lookup[, "Module"],
               stringsAsFactors=FALSE)
  for (curr_idx in 1:length(jsons_vec)) {
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

