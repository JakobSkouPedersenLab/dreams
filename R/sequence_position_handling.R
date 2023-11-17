#' Adjust Position Indices in BAM Data with CIGAR Strings
#'
#' @description This function takes a dataframe containing BAM file data and adjusts the position indices
#' based on CIGAR strings. CIGAR strings in BAM files provide information about alignment of reads to the
#' reference genome, including insertions, deletions, and padding. This function processes the CIGAR strings
#' to correct the read positions accordingly, accounting for insertions (I), deletions (D), hard clips (H),
#' and matches/mismatches (M). The adjusted position indices help in accurately mapping the reads back to
#' the reference sequence.
#'
#' @param df A dataframe containing BAM file data, which must include a column for CIGAR strings and a column
#' for position indices (`pos_idx`).
#'
#' @return A modified dataframe with corrected position indices accounting for CIGAR string operations.
#' The dataframe also includes additional columns that detail the size and index of insertions and deletions,
#' along with a logical column indicating whether a position is within a deletion.
#'
#' @keywords internal
#'
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#' @importFrom purrr map pmap pmap_dbl pmap_lgl
#' @importFrom stringr str_extract_all str_remove_all str_locate_all
#'
correct_pos_idx_w_cigar <- function(df) {
  df %>% mutate(

    # Extract numerical values from CIGAR strings and convert them to numeric
    cigar_pos = str_extract_all(string = .data$cigar, pattern = "\\d+(?=[HMID]+)") %>% lapply(as.numeric),

    # Identify indices of insertions in CIGAR strings
    insert_idx = str_remove_all(string = .data$cigar, pattern = "[\\d]") %>%
      str_locate_all(pattern = "I") %>%
      map(function(x) x[, "start"]),

    # Identify indices of deletions and hard clippings in CIGAR strings
    del_idx = .data$cigar %>%
      str_remove_all(pattern = "[\\d]") %>%
      str_locate_all(pattern = "[HD]") %>%
      map(function(x) x[, "start"]),

    # Calculate sizes of insertions
    ind_sz = map2(.data$cigar_pos, .data$insert_idx, function(x, y) x[y]),

    # Calculate sizes of deletions and hard clippings
    del_sz = map2(.data$cigar_pos, .data$del_idx, function(x, y) x[y]),

    # Compute the read index for insertions
    insert_read_index =
      map2(
        .data$cigar_pos, .data$insert_idx,
        function(cigar_pos, insert_idx) {
          cigar_pos[insert_idx] <- 0
          cumsum(cigar_pos)[insert_idx]
        }
      ),

    # Compute the read index for deletions and hard clippings
    del_read_index = pmap(
      list(
        .data$cigar_pos,
        .data$del_idx,
        .data$insert_idx
      ),
      function(cigar_pos, del_idx, insert_idx) {
        cigar_pos[insert_idx] <- 0
        cumsum(cigar_pos)[del_idx]
      }
    ),
    # Adjust position index for insertions
    ind_cor = pmap_dbl(
      list(
        .data$ind_sz,
        .data$insert_read_index,
        .data$pos_idx
      ),
      function(ind_sz, insert_read_index, pos_idx) {
        sum(ind_sz[insert_read_index < pos_idx])
      }
    ),

    # Adjust position index for deletions and hard clippings
    del_cor = pmap_dbl(
      list(
        .data$del_sz,
        .data$del_read_index,
        .data$pos_idx
      ),
      function(del_sz, del_read_index, pos_idx) {
        sum(del_sz[del_read_index < pos_idx])
      }
    ),

    # Determine if the position is within a deletion range
    is_in_deletion = pmap_lgl(
      list(
        .data$del_sz,
        .data$del_read_index,
        .data$pos_idx
      ),
      function(del_sz, del_read_index, pos_idx) {
        any((del_read_index - del_sz + 1 <= pos_idx) & (pos_idx <= del_read_index))
      }
    ),

    # Preserve original position index for reference
    old_pos_idx = .data$pos_idx,

    # Calculate the final adjusted position index
    pos_idx = .data$pos_idx + .data$ind_cor - .data$del_cor
  )
}


#' Get list of genomic positions with mismatches
#'
#' @description This function computes the genomic positions of mismatches in
#'   sequencing reads based on the MD tag from BAM files. The MD tag describes
#'   how the read aligns to the reference genome, including mismatches. The
#'   function identifies these mismatches and computes their positions relative
#'   to the start of the alignment, aiding in variant calling and analysis of
#'   sequencing errors.
#'
#' @param pos Genomic position
#' @param MDtag MD tag string from the BAM file that contains information about
#'   matches, mismatches, and deletions in the read alignment.
#'
#' @return A list of genomic positions where mismatches occur, corresponding to
#'   each MD tag provided.
#' @keywords internal
#' @importFrom stringr str_extract_all
#' @importFrom stringr str_count
#' @importFrom stringr str_which
#' @importFrom purrr map2
#' @importFrom rlang .data
get_mismatch_genomic_pos_list <- function(pos, MDtag) {
  # Extract the number of bases matching before each mismatch or
  # deletion event, convert to numeric
  leading_matches <- str_extract_all(MDtag, "\\d+(?=[\\^ATCGN]+)") %>% lapply(as.numeric)

  # Count the number of mismatched bases (or bases in a deletion event)
  # following each count of matching bases
  event_length <- str_extract_all(MDtag, "[\\^ATCGN]+") %>% lapply(str_count, pattern = "[ATCGN]")

  # Generate a cumulative position index list for the start of each event (match
  # or mismatch)
  pos_idx_list <- map2(leading_matches, event_length, function(x, y) cumsum(x + y))

  mismatch_idx <-
    str_extract_all(MDtag, "[\\^ATCGN]+") %>%
    lapply(str_which, pattern = "\\^", negate = TRUE)

  # Calculate the genomic positions of mismatches by adding the mismatch indices
  # to the start position The `- 1` corrects the offset to match the genomic
  # coordinate system
  mismatch_genomic_positions_list <-
    map2(pos_idx_list, mismatch_idx, function(x, y) x[y]) %>%
    map2(pos, function(x, y) x + y - 1)
  return(mismatch_genomic_positions_list)
}

#' Get_match_genomic_pos_list
#'
#' @param pos position_list
#' @param cigar cigar string list
#' @param MDtag mdtag list
#'
#' @return match_positions
#' @keywords internal
get_match_genomic_pos_list <- function(pos, cigar, MDtag) {
  # Remove insert from cigar
  cigar_inserts_removed <- cigar %>%
    str_remove_all(pattern = "[\\d]+I") %>%
    str_remove_all(pattern = "^[0-9]+[H]") %>%
    str_remove_all(pattern = "[0-9]+[H]$")

  event_lengths <-
    cigar_inserts_removed %>%
    str_extract_all(pattern = "\\d+(?=[MD]+)") %>%
    lapply(as.numeric)

  genomic_offset <- map(event_lengths, function(x) 1:sum(x))

  del_idx <-
    str_extract_all(cigar_inserts_removed, "[MD]") %>%
    lapply(str_which, pattern = "D")

  del_start <-
    map(event_lengths, cumsum) %>%
    map2(event_lengths, function(x, y) x - y + 1) %>%
    map2(del_idx, function(x, y) x[y])

  del_end <-
    map(event_lengths, cumsum) %>%
    map2(del_idx, function(x, y) x[y])

  del_positions <- map2(del_start, del_end, function(x, y) mapply(":", x, y) %>% unlist())
  mismatch_genomic_pos <- get_mismatch_genomic_pos_list(pos = pos, MDtag = MDtag)

  match_positions <-
    genomic_offset %>%
    # Remove deletions
    map2(del_positions, function(x, y) setdiff(x, y)) %>%
    # Get genomic position by adding read start position
    map2(pos, function(x, y) x + y - 1) %>%
    # Filter mismatch positions
    map2(mismatch_genomic_pos, function(x, y) setdiff(x, y))

  return(match_positions)
}


#' Sample_negative_read_positions
#'
#' @param bam_df bam dataframe from load_BAM
#' @param n_samples number of samples
#'
#' @return negative_read_positions
#' @keywords internal
#'
#' @importFrom purrr map_dbl
sample_negative_read_positions <- function(bam_df, n_samples) {

  # Sample random reads according to fragment size
  random_read_positions_df <-
    bam_df %>%
    sample_n(size = n_samples, weight = abs(.data$isize), replace = T)

  # Sample random negative position per read
  negative_read_positions <-
    random_read_positions_df %>%
    mutate(
      match_genomic_pos_list = get_match_genomic_pos_list(pos = .data$pos, cigar = .data$cigar, MDtag = .data$MD),
      genomic_pos = map_dbl(.data$match_genomic_pos_list, function(x) sample(x = x, size = 1))
    ) %>%
    select(-.data$match_genomic_pos_list)

  return(negative_read_positions)
}

#' Extract mismatch positions from BAM data
#'
#' @description This function filters a dataframe for reads that contain
#'   mismatches when compared to the reference genome, as indicated by the MD
#'   tag within the BAM data. It then extracts the genomic positions of these
#'   mismatches, providing a detailed account of where the sequencing reads do
#'   not match the reference sequence.
#'
#' @param bam_df A dataframe obtained from the `load_BAM` function.
#'
#' @return A dataframe of mismatch positions, which includes the genomic
#'   positions of all identified mismatches across the reads.
#' @keywords internal
#' @importFrom tidyr unnest
#' @importFrom dplyr filter
extract_mismatch_positions <- function(bam_df) {
  mismatch_positions <-
    bam_df %>%
    # Get all reads with at least one mismatch
    filter(str_detect(.data$MD, "\\d+[ATCG]+")) %>%
    # Select reads that have mismatches, which are indicated by the presence of
    # base letters (ATCG) in the MD tag
    mutate(
      genomic_pos = get_mismatch_genomic_pos_list(pos = .data$pos, MDtag = .data$MD)
    ) %>%
    # Expand the list column so that each mismatch position has its own row in
    # the dataframe.
    unnest(.data$genomic_pos)

  return(mismatch_positions)
}











