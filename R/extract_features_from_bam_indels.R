

#' Extract Features from BAM Data
#'
#' @param bam_df A DataFrame originating from 'load_BAM' and processed by extract_indel_info' containing alignment data.
#' @param reference_path Path to the reference genome file in FASTA format.
#' @param add_umi_features Check if umi information is available.
#'
#' @return A DataFrame with extracted features, including read positions and possibly UMI data.
#' @keywords internal
#'
#' @importFrom purrr map2_int
#'
#'
extract_features_from_bam_indels <- function(bam_df, reference_path, add_umi_features = all(c("cd", "ce") %in% colnames(bam_df))) {

  # If UMI features are asked for but not present
  if (add_umi_features & !all(c("cd", "ce") %in% colnames(bam_df))) {
    stop("UMI features (ce and cd) are not available in bam_df!")
  }
  # Return an empty dataframe if the input dataframe has no rows
  if (nrow(bam_df) == 0) {
    return(data.frame())
  }

  # Create genomic position features by selecting distinct chromosomes, positions, and strands
  # and calculating various genomic features such as local complexity and GC content
  genomic_pos_feature_df <-
    distinct(
      base::data.frame(
        chr = bam_df$chr,
        genomic_pos = bam_df$genomic_pos,
        strand = bam_df$strand
      )
    ) %>%
    mutate(
      context11 =
        get_reference_seq(
        chr = .data$chr,
        genomic_pos = .data$genomic_pos,
        buffer = 5,
        reference_path = reference_path
      ),
      local_complexity_1 = calc_string_entropy_k_mer(.data$context11, k = 1),
      local_complexity_2 = calc_string_entropy_k_mer(.data$context11, k = 2),
      local_GC = str_count(.data$context11, "[GC]") / 11,
      ctx_minus1 = substring(.data$context11, 5, 5),
      ref = substring(.data$context11, 6, 6),
      ctx_plus1 = substring(.data$context11, 7, 7),
      trinucleotide_ctx = paste0(.data$ctx_minus1, .data$ref, .data$ctx_plus1)
    )

  # Create read-specific features from the BAM data, adjusting for insertions, deletions, and errors
  read_feature_df <-
    bam_df %>%
    mutate(
      pos_idx = .data$genomic_pos - .data$pos + 1,
      corrected_seq = mapply(correct_seq, .data$cigar, .data$seq),
      fragment_size = abs(.data$isize),
      seq_length = nchar(.data$seq),
      read_index = if_else(.data$strand == "fwd", .data$pos_idx, nchar(.data$corrected_seq) - .data$pos_idx + 1),
      first_in_pair = as.integer(as.logical(bitwAnd(.data$flag, 64))),
      n_errors_in_read = str_count(.data$MD, "\\d+[ATCG]"),
      n_insertions_in_read = str_count(.data$cigar, "I"),
      n_deletions_in_read = str_count(.data$cigar, "D")
    ) %>%
    # TODO: Move to filter function! Or do before calling this function!
    filter(.data$fragment_size != 0)

  # Define the features to be included in the output
  selected_features <-
    c(
      "qname", "chr", "genomic_pos", "obs", "ref",
      "strand", "first_in_pair", "read_index", "fragment_size",
      "ctx_minus1", "ctx_plus1", "trinucleotide_ctx", "context11",
      "local_complexity_1", "local_complexity_2", "local_GC",
      "n_insertions_in_read", "n_deletions_in_read", "seq_length"
    )


  # Add UMI features if asked
  if (add_umi_features) {
    read_feature_df <-
      read_feature_df %>%
      mutate(
        umi_count = map2_int(.data$cd, .data$pos_idx, function(cd, pos_idx) cd[pos_idx]),
        umi_errors = map2_int(.data$ce, .data$pos_idx, function(ce, pos_idx) ce[pos_idx])
      )

    # Add UMI features to selection
    selected_features <- c(selected_features, c("umi_count", "umi_errors"))
  }

  # Join and select features: Read, genomic positions and UMI
  feature_df <-
    left_join(read_feature_df, genomic_pos_feature_df,
      by = c("chr", "genomic_pos", "strand")
    ) %>%
    mutate(
      obs = case_when(
        substring(corrected_seq, pos_idx, pos_idx) == "D" ~ "D",
        substring(corrected_seq, pos_idx, pos_idx) == "I" ~ "I",
        substring(corrected_seq, pos_idx, pos_idx) %in% c("A", "C", "T", "G") ~ substring(corrected_seq, pos_idx, pos_idx),
        TRUE ~ "N"
      )) %>%
    select(all_of(selected_features))
  return(feature_df)
}



#' Correct a DNA Sequence Based on a CIGAR String
#'
#' This function corrects a DNA sequence based on a provided CIGAR string. It first expands the CIGAR string,
#' then inserts 'D's into the DNA sequence at positions corresponding to 'D's in the expanded CIGAR string.
#' It also replaces characters in the DNA sequence with 'I's at positions corresponding to 'I's in the CIGAR string.
#' After applying these modifications, the function calls `clean_insertions` to clean up the sequence.
#'
#' @param cigar A string representing the CIGAR sequence. Each character in the CIGAR string (such as M, I, D)
#' represents a type of alignment operation.
#' @param sequence A string representing the DNA sequence that needs to be corrected based on the CIGAR string.
#'
#' @return A DNA sequence string that has been corrected based on the CIGAR string. The function returns the
#' sequence after applying deletions, insertions, and further cleaning through the `clean_insertions` function.
#'
correct_seq <- function(cigar, sequence) {

  #Expand cigar
  expanded_cigar <- expand_cigar(cigar)

  # Find the positions of 'D' in the first string
  positions <- which(strsplit(expanded_cigar, "")[[1]] == "D")

  # Insert 'D's into the second string at the corresponding positions
  for (pos in positions) {
    sequence <- paste0(substr(sequence, 1, pos - 1), "D", substr(sequence, pos, nchar(sequence)))
  }
  # Convert both the CIGAR string and the sequence into character vectors
  cigar_chars <- strsplit(expanded_cigar, "")[[1]]
  sequence_chars <- strsplit(sequence, "")[[1]]

  # Find the positions of 'I' in the CIGAR string
  I_positions <- which(cigar_chars == "I")

  # Replace corresponding positions in the sequence with 'I'
  sequence_chars[I_positions] <- "I"

  # Combine the modified sequence into a single string and clean insertions
  corrected_seq <- clean_insertions(paste(sequence_chars, collapse = ""))

  return(corrected_seq)
}

