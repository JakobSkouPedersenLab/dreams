#' Expand CIGAR string
#'
#' @param cigar CIGAR string
#'
#' @return List of expanded CIGAR string
#' @keywords internal
#'

expand_cigar <- function(cigar) {

  cigar = str_remove_all(cigar, pattern = "^[0-9]+[HD]")
  cigar = str_remove_all(cigar, pattern = "[0-9]+[HD]$")
  # Extract the operation counts and types
  operations <- str_extract_all(cigar, "[MIDNSHP=X]")[[1]]
  counts <- as.numeric(str_extract_all(cigar, "\\d+")[[1]])

  # Use sapply to repeat each operation count times and concatenate them
  expanded_cigar <- paste0(sapply(seq_along(operations), function(i) {
    strrep(operations[i], counts[i])
  }), collapse = "")

  return(expanded_cigar)
}

#' Expand MD tag
#'
#' @param md MD tag (`string`)
#'
#' @return List of expanded MD tag
#' @keywords internal
#'
expand_md <- function(md) {
  # Use strsplit with a regular expression that matches numbers or single non-digit characters
  split_string <- strsplit(md, "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)", perl = TRUE)
  split_elements <- split_string[[1]]

  # Use a vectorized approach to create a string with 'M' for numbers and 'D' for '^'
  str <- vapply(split_elements, function(element) {
    if (grepl("^[0-9]+$", element)) {
      return(strrep("M", as.numeric(element)))
    } else if (startsWith(element, "^")) {
      return(strrep("D", nchar(element) - 1))
    } else {
      return(element)
    }
  }, character(1))

  # Collapse the vector into a single string
  str <- paste0(str, collapse = "")
  return(str)
}

#' Combined CIGAR and MD tag
#'
#' @param cigar CIGAR string
#' @param md MD tag (`string`)
#'
#' @return List of combined CIGAR and MD tag
#' @keywords internal
#'
combine_sequences <- function(cigar, md) {
  cigar <- expand_cigar(cigar)
  md <- expand_md(md)
  # Split both strings into arrays of characters
  cigar_chars <- strsplit(cigar, "")[[1]]
  md_chars <- strsplit(md, "")[[1]]

  # Use vectorized comparison to identify positions where md_chars is not 'M'
  # and it should override the cigar_chars
  replace_indices <- which(md_chars != "M")

  # Create a combined vector initialized to cigar_chars
  combined_chars <- cigar_chars

  # Replace the corresponding indices in combined_chars with md_chars
  combined_chars[replace_indices] <- md_chars[replace_indices]

  # Collapse the character vector back into a single string
  combined <- paste0(combined_chars, collapse = "")
  return(combined)
}


