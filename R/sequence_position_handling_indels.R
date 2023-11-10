
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
  # Find all matches for 'I's.
  I_positions <- gregexpr("I", cigar)
  # Convert the positions to a list of more readable numbers (1-based indices).
  I_positions <- unlist(I_positions)
  # Remove the -1 values which indicate no match found.
  I_positions <- I_positions[I_positions != -1]

  # Initialize an empty vector to store the indexes for D's.
  D_positions <- c()
  # Loop through the string
  for (i in seq(nchar(cigar))) {
    # Check if the character is 'D'
    if (substr(cigar, i, i) == "D") {
      # Check if it's the start of a 'D' sequence or if the previous character
      # is not a 'D'
      if (i == 1 || substr(cigar, i - 1, i - 1) != "D") {
        D_positions <- c(D_positions, i)
      }
    }
  }
  # Return a list containing positions of 'I's and the specified 'D's. #
  # Calculate the genomic positions of mismatches by adding the mismatch indices
  # to the start position The `- 1` corrects the offset to match the genomic
  # coordinate system
  list(I_positions = I_positions + pos-1, D_positions = D_positions + pos-1)
}



#' Combine Positions of Insertions and Deletions
#'
#' @description This function identifies and combines the positions of
#'   insertions ('I') and deletions ('D'). The positions are returned as a
#'   single vector, calculated based on the start position of the CIGAR string
#'   in the genome.
#'
#' @param pos Integer, the start position of the CIGAR string in the genome.
#' @param cigar String, the CIGAR string representing genomic alignments.
#'
#' @return A numeric vector containing the combined positions of insertions and
#'   deletions. Returns an empty vector if there are no insertions or deletions.
#' @keywords internal
#' @examples
#' # Assuming you have a pos and a cigar string
#' # pos <- 10
#' # cigar <- "8M1I4M2D"
#' # combined_positions <- combine_positions(pos, cigar)
#'
combine_positions <- function(pos, cigar) {
  # Obtain the positions of insertions and deletions using the
  # get_positions_of_indels function
  positions <- get_positions_indels(pos, cigar)

  # Combine the positions of 'I's and the first 'D' in each sequence of 'D's
  # into a single vector and return this combined list of positions
  c(positions$I_positions, positions$D_positions)
}


#' Extract lengths of insertions and deletions
#'
#' This function extracts numeric lengths that appear before 'I' and 'D' in a
#' CIGAR string.
#'
#' @param cigar A character string representing a CIGAR string.
#' @return A list containing two elements: `I_length` and `D_length`.
#' Each element is a vector of numbers (as character strings) that appear before 'I' and 'D' in the CIGAR string.
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

extract_lengths_insertions <- function(cigar){
  return(extract_lengths_indels(cigar)$I_length)
}

extract_lengths_deletions <- function(cigar){
  return(extract_lengths_indels(cigar)$D_length)
}




#' Extract indels positions from BAM data
#'
#' @description This function filters a dataframe for reads that contain indels
#'   when compared to the reference genome, as indicated by the cigar string in
#'   the BAM data. It then extracts the genomic positions of these indels,
#'   providing a detailed account of where the sequencing reads do not match the
#'   reference sequence.
#'
#' @param bam_df A dataframe obtained from the `load_BAM` function.
#'
#' @return A dataframe of indel positions, which includes the genomic
#'   positions of all identified indels across the reads.
#' @keywords internal
#' @importFrom tidyr unnest
extract_mismatch_positions_indels <- function(bam_df) {
  indels_positions <-
    bam_df %>%
    # Get all reads with at least one indel
    filter(str_detect(.data$cigar, "I|D"))
  # len_insertions <- unlist(lapply(indels_positions$cigar, extract_lengths_insertions))
  indels_positions <- indels_positions %>%
    # Select reads that have indels and provide the lengths.
    mutate(genomic_pos = combine_positions(pos = .data$pos, cigar = .data$cigar)) %>%
    # Expand the list column so that each mismatch position has its own row in
    # the dataframe.
    unnest(genomic_pos)


  return(indels_positions)
}




