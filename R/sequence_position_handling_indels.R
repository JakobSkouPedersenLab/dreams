
#' Expand CIGAR String into Sequence of Operations
#'
#' @description This function expands a condensed CIGAR string into a sequence
#'   of operations.
#' @param cigar A character string representing the CIGAR operations from an
#'   alignment.
#'
#' @return A character string with each CIGAR operation expanded to show the
#'   full sequence of alignment operations. For example, '2M1I' would expand to
#'   'MMI'.
#' @keywords internal
#'
expand_cigar <- function(cigar) {

  # Extract the operation counts and types
  operations <- str_extract_all(cigar, "[MIDNP=X]")[[1]]
  counts <- as.numeric(str_extract_all(cigar, "\\d+")[[1]])

  # Use sapply to repeat each operation count times and concatenate them
  expanded_cigar <- paste0(sapply(seq_along(operations), function(i) {
    strrep(operations[i], counts[i])
  }), collapse = "")

  return(expanded_cigar)
}

#' Clean Insertions in CIGAR Sequence
#'
#' @description This function adjusts the sequence of CIGAR operations by
#'   cleaning up insertion operations. It removes any insertions that follow
#'   deletions and consolidates consecutive insertions that follow any operation
#'   other than deletions into a single insertion, while removing the operation
#'   that precedes the insertions.
#' @param expanded_cigar A string representing the sequence output from
#'   expand_cigar.
#'
#' @return A string of the CIGAR sequence after cleaning up the insertion
#'   operations. For example, 'MMDDIMMII' would be to MMDDMI'.
#' @keywords internal
#'
clean_insertions <- function(expanded_cigar) {

   # Remove insertions that follow deletions.
  cigar <- gsub("DI+", "D", expanded_cigar)

  # For 'I's following any letter except 'D', keep one 'I' and remove the
  # preceding letter
  cleaned_cigar <- gsub("([A-HJ-Z])I+", "I", cigar)

  return(cleaned_cigar)
}

#' Convert Expanded CIGAR String to Standard CIGAR Format
#'
#' @description This function takes an expanded CIGAR string (where each
#'   operation is represented individually, e.g., "MMMMMMMDDDDDDDMIDD") and
#'   converts it into standard CIGAR format (where consecutive operations are
#'   counted, e.g., "7M7D1M1I2D").
#'
#' @param cleaned_cigar A character string representing the expanded CIGAR
#'   string.
#'
#' @return A character string in standard CIGAR format.
#'
#' @keywords internal
#'
convert_to_cigar <- function(cleaned_cigar) {
  # Use regular expressions to find matches of consecutive characters
  matches <- gregexpr("(.)\\1*", cleaned_cigar, perl = TRUE)[[1]]

  # Get the length of each match
  lengths <- attr(matches, "match.length")

  # Extract the character for each match
  characters <- sapply(matches, function(i) substr(cleaned_cigar, i, i))

  # Combine lengths and characters to form the standard CIGAR string
  cigar <- paste0(lengths, characters, collapse = "")

  return(cigar)
}

#' Extract Start Positions of Indels from a CIGAR String
#'
#' @description This function analyzes a CIGAR string and extracts the start
#'   positions of insertions and deletions. It parses the CIGAR string,
#'   calculates the cumulative lengths of operations, and then identifies the
#'   start positions of the 'I' and 'D' segments.
#'
#' @param cigar A character string representing the CIGAR format.
#'
#' @return An integer vector containing the start positions of each insertion
#'   and deletion in the CIGAR string. Positions are 1-based.
#'
get_indels_start_positions <- function(cigar){

  # Extract all segments of the CIGAR string
  segments <- regmatches(cigar, gregexpr("\\d+[MDI]", cigar))[[1]]

  # Separate numbers and characters, and compute cumulative lengths
  numbers <- as.numeric(sub("[MDI]", "", segments))
  chars <- gsub("[0-9]", "", segments)
  cum_lengths <- cumsum(numbers)

  # Identify segments that are I or D and calculate start positions
  id_positions <- cum_lengths[chars %in% c("I", "D")] - numbers[chars %in% c("I", "D")] + 1

  return(id_positions)
}



#' Extract Indel Information from CIGAR String and Genomic Position
#'
#' This function takes a genomic position and a CIGAR string as input and
#' extracts information about indels. It processes the CIGAR string to determine
#' the genomic positions, lengths, and indel type.
#'
#' @param pos An integer representing the starting genomic position of the CIGAR
#'   string.
#' @param cigar A character string representing the CIGAR format.
#'
#' @return A list containing the following elements:
#'         - `genomic_pos`: An integer vector of the genomic positions of the indels.
#'         - `indel_length`: An integer vector of the lengths of the indels.
#'         - `indel_type`: A character vector of the types of the indels
#'                        (either 'I' for insertion or 'D' for deletion).
#'
#'
get_indel_info <- function(pos, cigar) {

  # Expand the CIGAR string to include full representation of each operation
  expanded_cigar <- expand_cigar(cigar)

  # Clean the expanded CIGAR
  cleaned_cigar <- clean_insertions(expanded_cigar)

  convert_to_cigar <- convert_to_cigar(cleaned_cigar)

  # Extract genomic positions of indels
  genomic_pos <- get_indels_start_positions(convert_to_cigar)

  # Extract matches for indel length and type from the original CIGAR string
  matches <- regmatches(cigar, gregexpr("\\d+[ID]", convert_to_cigar, perl = TRUE))[[1]]

  # Calculate the length of each indel
  indel_length <- as.numeric(sub("[ID]", "", matches))

  # Return a list with genomic positions, lengths, and types of indels
  list(genomic_pos = genomic_pos + pos-1,
       indel_length = indel_length)
}


#' Extract Genomic Positions of Indels from a CIGAR String
#'
#' This function uses the `get_indel_info` function to extract the genomic
#' positions of indels from a given CIGAR string, starting from a specified
#' genomic position.
#'
#' @param pos Integer, the starting genomic position for the CIGAR string
#'   analysis.
#' @param cigar String, the CIGAR string representing genomic alignments.
#'
#' @return An integer vector containing the genomic positions of the indels in
#'   the CIGAR string.
#'
#' @keywords internal
#'
get_indel_genomic_pos <- function(pos, cigar){
  return(get_indel_info(pos, cigar)$genomic_pos)
}

#' Extract Lengths of Indels from a CIGAR String
#'
#' This function uses the `get_indel_info` function to extract the lengths of
#' indels from a given CIGAR string, starting from a specified genomic position.
#'
#' @param pos Integer, the starting genomic position for the CIGAR string
#'   analysis.
#' @param cigar String, the CIGAR string representing genomic alignments.
#'
#' @return An integer vector containing the lengths of the indels in the CIGAR
#'   string.
#'
#' @keywords internal
#'
get_indel_length <- function(pos, cigar){
  return(get_indel_info(pos, cigar)$indel_length)
}

#' Get Indel Length with Zero Default
#'
#' This function calculates the indel length based on the provided position (pos)
#' and cigar string (cigar). It wraps around the `get_indel_length` function and ensures
#' that a 0 is returned whenever the indel length is either `NULL` or 0.
#'
#' @param pos A numerical value or vector representing the position(s).
#' @param cigar A string or vector of strings representing the CIGAR data.
#'
#' @return A numerical value or vector representing the indel lengths.
#'         If the indel length calculated by `get_indel_length` is `NULL` or 0,
#'         this function returns 0 instead.
#'
get_indel_length_with_zero <- function(pos, cigar) {
  indel_length_result <- get_indel_length(pos, cigar)
  if (identical(indel_length_result, numeric(0))) {
    indel_length_result = 0
  }
  return(indel_length_result)
}


#' Extract Indel Information from BAM Data Frame
#'
#' This function processes a data frame containing BAM file information to
#' extract indel information. It filters rows with insertions or deletions,
#' calculates indel lengths and types, and generates genomic positions and indel
#' sequences.
#'
#' @param bam_df A data frame containing BAM file information, including columns
#'   for position (`pos`), CIGAR string (`cigar`), and sequence (`seq`).
#'
#' @return A modified version of the input data frame `bam_df`, which includes
#'   additional columns for genomic position (`genomic_pos`), indel length
#'   (`indel_length`), indel type (`indel_type`), and indel sequence
#'   (`indel_seq`). Each row corresponds to an indel event.
#'
#' @importFrom dplyr filter mutate
#' @importFrom stringr str_detect
#' @importFrom tidyr unnest
#'
#' @keywords internal
#'
extract_indel_info <- function(bam_df) {

  indels <- bam_df %>%
    filter(str_detect(.data$cigar, "\\d+[ID]"))

  if (nrow(indels) == 0) {
    return(data.frame())
  }

  indels <- indels %>%
    mutate(genomic_pos = Map(get_indel_genomic_pos, .data$pos, .data$cigar)) %>%
    unnest(genomic_pos)


  return(indels)
}





