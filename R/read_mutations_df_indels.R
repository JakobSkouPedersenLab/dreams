#' Read and Modify Mutations Data Frame
#'
#' This function processes a data frame containing mutation information. It modifies the 'ref', 'alt', and 'genomic_pos' columns based on certain conditions.
#' Specifically, if the 'ref' field contains more than one character, the function retains only the second character, marks the mutation as a deletion (by setting 'alt' to "D"),
#' and increments the 'genomic_pos' by 1. If the 'alt' field contains more than one character, the mutation is marked as an insertion (by setting 'alt' to "I").
#'
#' @param mutations_df A data frame representing mutations. Each row should contain the fields 'ref', 'alt', and 'genomic_pos'.
#'
#' @return A modified data frame with updated mutation information.
#'
#' @keywords internal
#'
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
