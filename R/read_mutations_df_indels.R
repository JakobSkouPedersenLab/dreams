read_mutations_df <- function(mutations_df){
  for (i in 1:nrow(mutations_df)){
    if (nchar(mutations_df[i, "ref"]) > 1){
      # Keep only the second character in ref
      mutations_df[i, "ref"] = substring(mutations_df[i, "ref"], 2, 2)
      # Mark as deletion
      mutations_df[i, "alt"] = "D"
      mutations_df[i, "genomic_pos"] = mutations_df[i, "genomic_pos"] + 1
    } else if (nchar(mutations_df[i, "alt"]) > 1)
      # Mark as insertion
      mutations_df[i, "alt"] = "I"
  }
  return(mutations_df)
}
