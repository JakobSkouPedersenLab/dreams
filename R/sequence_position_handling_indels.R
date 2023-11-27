
#' Expand CIGAR String into Sequence of Operations
#'
#' @description This function expands a condensed CIGAR string into a sequence
#'   of operations. It removes operations at the beginning and end of the string
#'   that represent hard clipping ('H') and deletions ('D'), as these do not
#'   contribute to the alignment sequence.
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

#' Get genomic positions of insertions and deletions
#'
#' @description This function takes a CIGAR string as input and returns the
#'   genomic start positions of insertions and deletions. Deletions start on
#'   their actual position, whereas insertions start positions are annotated on
#'   the base prior to the insertion.
#' @param pos Start position of CIGAR string in genome.
#' @param cigar CIGAR string.
#'
#' @return A list with two elements:
#'         - `I_positions`: a numeric vector with the positions of 'I's.
#'         - `D_positions`: a numeric vector with the positions of the first 'D'
#'   in each sequence of 'D's.
#'
#' @keywords internal
get_positions_indels <- function(pos, cigar) {
  expanded_cigar <- expand_cigar(cigar)
  cigar <- clean_insertions(expanded_cigar)

  # Find all start positions of insertions
  I_positions <- as.numeric(unlist(gregexpr("I+", cigar, perl = TRUE)))
  # Filter out all values smaller than 1
  I_positions <- I_positions[I_positions >= 1]

  # Find start positions of deletions
  D_positions <- as.numeric(unlist(gregexpr("D+", cigar, perl = TRUE)))
  # Filter out all values smaller than 1
  D_positions <- D_positions[D_positions >= 1]
  # Return a list containing positions of 'I's and 'D's.
  list(I_positions = I_positions + pos-1, D_positions = D_positions + pos-1)
}

#' Combine Positions from CIGAR Strings
#'
#' This function takes two vectors `pos` and `cigar`, and combines the positions
#' of insertions ('I') and the first deletion ('D') in each sequence of 'D's from
#' the CIGAR strings. It returns a list of combined positions for each pair of
#' `pos` and `cigar`. The function relies on `get_positions_indels` to extract
#' insertion and deletion positions from CIGAR strings.
#'
#' @param pos A numeric vector representing positions.
#' @param cigar A character vector of CIGAR strings. Each CIGAR string
#'   corresponds to a position in `pos`. The length of `cigar` must match the
#'   length of `pos`.
#'
#' @return A list where each element is a numeric vector. Each vector contains
#'   combined positions of 'I's and the first 'D' from each sequence of 'D's in
#'   the corresponding CIGAR string.
#'
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#'
#' @keywords internal
#'
combine_positions <- function(pos, cigar) {

  combined_positions_list <- lapply(seq_along(pos), function(i) {
    pos <- pos[i]
    cigar <- cigar[i]

    positions <- get_positions_indels(pos, cigar)


    c(positions$I_positions, positions$D_positions)
  })
  return(combined_positions_list)
}


#' Extract lengths of insertions and deletions
#'
#' This function extracts numeric lengths that appear before 'I' and 'D' in a
#' CIGAR string.
#'
#' @param cigar A character string representing a CIGAR string.
#' @return A list containing two elements: `I_length` and `D_length`. Each
#'   element is a vector of numbers (as character strings) that appear before
#'   'I' and 'D' in the CIGAR string.
#' @keywords internal
#' @importFrom stringr str_extract_all
#'
extract_lengths_indels <- function(cigar) {
  # Use stringr to extract numbers before 'I'
  numbers_D_0 <- gsub("(\\d+)(?=D)", "0", cigar, perl = TRUE)
  I_length <- as.numeric(str_extract_all(numbers_D_0,  "\\d+(?=[ID])")[[1]])

  # Use stringr to extract numbers before 'D'
  numbers_I_0 <- gsub("(\\d+)(?=I)", "0", cigar, perl = TRUE)
  D_length <- as.numeric(str_extract_all(numbers_I_0,  "\\d+(?=[ID])")[[1]])

  # Return the results as a list
  list(I_length = I_length, D_length = D_length)
}

#' Extract Lengths of Insertions from CIGAR Strings
#'
#' @description
#' This function processes a given CIGAR string and extracts the lengths of all
#' insertion operations. It leverages the `extract_lengths_indels` function
#' to parse the CIGAR string and isolate the lengths corresponding to insertions.
#'
#' @param cigar A character vector of CIGAR strings, each representing the alignment
#'   of a single read in a BAM file.
#'
#' @return A numeric vector where each element represents the total length of
#'   insertions in the corresponding CIGAR string.
#' @keywords internal
#'
extract_lengths_insertions <- function(cigar) {
  # Extract lengths of insertions from the CIGAR string using
  # the extract_lengths_indels function and return the result.
  return(extract_lengths_indels(cigar)$I_length)
}

#' Extract Lengths of Deletions from CIGAR Strings
#'
#' @description
#' This function processes a given CIGAR string and extracts the lengths of all
#' deletion operations. It utilizes the `extract_lengths_indels` function
#' to parse the CIGAR string and isolate the lengths corresponding to deletions.
#'
#' @param cigar A character vector of CIGAR strings, each representing the alignment
#'   of a single read in a BAM file.
#'
#' @return A numeric vector where each element represents the total length of
#'   deletions in the corresponding CIGAR string.
#' @keywords internal
#' @seealso \code{\link{extract_lengths_indels}}
#'
extract_lengths_deletions <- function(cigar) {
  # Extract lengths of deletions from the CIGAR string using
  # the extract_lengths_indels function and return the result.
  return(extract_lengths_indels(cigar)$D_length)
}

#' Extract Indel Positions and lengths from BAM Data
#'
#' @description This function processes a dataframe obtained from `load_BAM`cto
#' identify and extract the genomic positions of indels. It filters the
#' dataframe for reads that contain indels, as indicated by the CIGAR string,
#' and then calculates the precise genomic positions, length and bases of these indels.
#'
#' @param bam_df A dataframe obtained from the `load_BAM` function, containing
#'   BAM data including fields for CIGAR strings and positions.
#'
#' @return A dataframe detailing the genomic positions of indels. This includes
#'   columns for the positions of indels and their respective
#'   lengths and base composition.
#'
#' @keywords internal
#'
#' @importFrom tidyr unnest
#' @importFrom dplyr filter
#' @importFrom dplyr mutate
#' @importFrom stringr str_detect
#'
extract_indel_pos_len <- function(bam_df) {
  # Filter to retain only reads with indels (based on CIGAR string)
  indels_positions <- bam_df %>%
    filter(str_detect(.data$cigar, "\\d+[ID]"))

  # Extract lengths of insertions and deletions from the CIGAR string
  len_insertions <- unlist(lapply(indels_positions$cigar, extract_lengths_insertions))
  len_deletions <- unlist(lapply(indels_positions$cigar, extract_lengths_deletions))

  # Add a column with genomic positions of indels
  indels_positions <- indels_positions %>%
    mutate(genomic_pos = combine_positions(pos = .data$pos, cigar = .data$cigar)) %>%
    # Expand each list of genomic positions into individual rows
    unnest(.data$genomic_pos)

  # Add columns for lengths of insertions and deletions
  indels_len <- indels_positions %>%
    mutate(len_insertions = len_insertions,
           len_deletions = len_deletions)

  indels_str <- indels_len %>%
    mutate(str_insertion = substring(.data$seq, .data$genomic_pos + 1 - .data$pos, .data$genomic_pos - .data$pos + .data$len_insertions),
           str_deletion = substring(.data$seq, .data$genomic_pos + 1 - .data$pos, .data$genomic_pos - .data$pos + .data$len_deletions))

  return(indels_str)
}





