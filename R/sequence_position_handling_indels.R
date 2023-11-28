
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
#' @param cigar A string representing the sequence output from expand_cigar.
#'
#' @return A string of the CIGAR sequence after cleaning up the insertion
#'   operations.
#' @keywords internal
#'
clean_insertions <- function(cigar) {

   # Remove insertions that follow deletions.
  cigar <- gsub("DI+", "D", cigar)

  # For 'I's following any letter except 'D', keep one 'I' and remove the
  # preceding letter
  cigar <- gsub("([A-HJ-Z])I+", "I", cigar)

  return(cigar)
}


#' Extract Indel Information from a CIGAR String
#'
#' This function takes a genomic position and a CIGAR string as input and extracts
#' information about indels present in the CIGAR string.
#' It returns the genomic positions, lengths, and types of these indels.
#'
#' @param pos Integer, the starting genomic position of the CIGAR string.
#' @param cigar String, the CIGAR string representing genomic alignments.
#'
#' @return A list containing three elements: `genomic_pos` ,
#'         `indel_length`, and `indel_type` .
#'
#' @importFrom stringr str_detect
#'
#' @keywords internal
#'
get_indel_info <- function(pos, cigar) {

  # Expand the CIGAR string to include full representation of each operation
  expanded_cigar <- expand_cigar(cigar)

  # Clean the expanded CIGAR
  cleaned_cigar <- clean_insertions(expanded_cigar)

  # Extract genomic positions of indels
  genomic_pos <- as.numeric(unlist(gregexpr("[ID]+", cleaned_cigar, perl = TRUE)))
  genomic_pos <- genomic_pos[genomic_pos >= 1]

  # Extract matches for indel length and type from the original CIGAR string
  matches <- regmatches(cigar, gregexpr("\\d+[ID]", cigar, perl = TRUE))[[1]]

  # Calculate the length of each indel
  indel_length <- as.numeric(sub("[ID]", "", matches))

  # Determine the type of each indel
  indel_type <- substring(matches, nchar(matches))

  # Return a list with genomic positions, lengths, and types of indels
  list(genomic_pos = genomic_pos + pos-1,
       indel_length = indel_length,
       indel_type = indel_type)
}


#' Extract Genomic Positions of Indels from a CIGAR String
#'
#' This function uses the `get_indel_info` function to extract the genomic positions of
#' indels from a given CIGAR string, starting from a specified genomic position.
#'
#' @param pos Integer, the starting genomic position for the CIGAR string analysis.
#' @param cigar String, the CIGAR string representing genomic alignments.
#'
#' @return An integer vector containing the genomic positions of the indels in the CIGAR string.
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
#' @param pos Integer, the starting genomic position for the CIGAR string analysis.
#' @param cigar String, the CIGAR string representing genomic alignments.
#'
#' @return An integer vector containing the lengths of the indels in the CIGAR string.
#'
#' @keywords internal
#'
get_indel_length <- function(pos, cigar){
  return(get_indel_info(pos, cigar)$indel_length)
}

#' Extract Types of Indels from a CIGAR String
#'
#' This function utilizes `get_indel_info` to determine the types
#' of indels in a given CIGAR string.
#'
#' @param pos Integer, the starting genomic position for the CIGAR string.
#' @param cigar String, the CIGAR string representing genomic alignments.
#'
#' @return A character vector containing the types of indels (either 'I' for insertions or 'D' for deletions)
#'         in the CIGAR string.
#'
#' @keywords internal
#'
get_indel_type <- function(pos, cigar){
  return(get_indel_info(pos, cigar)$indel_type)
}


#' Extract Indel Information from BAM Data Frame
#'
#' This function processes a data frame containing BAM file information to extract
#' indel information. It filters rows with insertions or deletions,
#' calculates indel lengths and types, and generates genomic positions and indel sequences.
#'
#' @param bam_df A data frame containing BAM file information, including columns for position (`pos`),
#'               CIGAR string (`cigar`), and sequence (`seq`).
#'
#' @return A modified version of the input data frame `bam_df`, which includes additional columns
#'         for genomic position (`genomic_pos`), indel length (`indel_length`), indel type (`indel_type`),
#'         and indel sequence (`indel_seq`). Each row corresponds to an indel event.
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

  indel_length = unlist(Map(get_indel_length, indels$pos, indels$cigar))
  indel_type = unlist(Map(get_indel_type, indels$pos, indels$cigar))

  indels <- indels %>%
    mutate(genomic_pos = Map(get_indel_genomic_pos, .data$pos, .data$cigar)) %>%
    unnest(genomic_pos)

  indels <- indels %>%
    mutate(
      indel_length = indel_length,
      indel_type = indel_type,
      indel_seq = substring(.data$seq, .data$genomic_pos + 1 - .data$pos, .data$genomic_pos - .data$pos + .data$indel_length)

    )

  return(indels)
}





