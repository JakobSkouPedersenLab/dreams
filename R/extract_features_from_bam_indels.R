

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
      pos_idx = .data$genomic_pos - .data$pos + 1
    ) %>%
    correct_pos_idx_w_cigar() %>%
    mutate(
      seq_corrected = mapply(remove_insertions, .data$cigar, .data$seq),
      cleaned_cigar = lapply(lapply(.data$cigar, expand_cigar), clean_insertions),
      fragment_size = abs(.data$isize),
      seq_length = nchar(.data$seq),
      read_index = if_else(.data$strand == "fwd", .data$pos_idx, - .data$pos_idx),
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
      "n_other_errors", "n_insertions_in_read", "n_deletions_in_read",
      "seq_length", "indel_length", "indel_type", "indel_seq"
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
        substring(cleaned_cigar, pos_idx, pos_idx) == "D" ~ "D",
        substring(cleaned_cigar, pos_idx, pos_idx) == "I" ~ "I",
        substring(cleaned_cigar, pos_idx, pos_idx) == "M" &
          substring(seq_corrected, pos_idx, pos_idx) == ref ~ substring(seq_corrected, pos_idx, pos_idx),
        TRUE ~ "N"
      ))
  feature_df <- feature_df %>%
    mutate(
      n_other_errors = .data$n_errors_in_read - ifelse((!.data$is_in_deletion) & (.data$obs != .data$ref), 1, 0)
    )%>%
    select(all_of(selected_features))
  return(feature_df)
}

#' Extract Features from BAM for Indels in Negative Samples
#'
#' This function processes a BAM file data frame and a reference path to extract features
#' specifically for indels (insertions or deletions) in negative samples. It enhances the data frame
#' with observed values set to reference values.
#'
#' @param bam_df A data frame representing BAM file reads.
#' @param reference_path Path to the reference used in feature extraction.
#' @return A modified data frame where each entry is enriched with additional features, particularly
#' for indels in negative samples. The 'obs' column is added, which is set equal to the 'ref' column.
#'
#' @keywords internal
#'
extract_features_from_bam_indels_negatives <- function(bam_df, reference_path){
  extract_features_negatives <- extract_features_from_bam_indels(
    bam_df,
    reference_path) %>%
    mutate(obs = ref)
  return(extract_features_negatives)
  }

#' Remove Insertions from Observed Sequence Based on CIGAR String
#'
#' This function processes an observed sequence based on a given CIGAR string,
#' removing any insertions (denoted as 'I' in the CIGAR string) from the observed sequence.
#' It first expands the CIGAR string to match the length of the observed sequence
#' and then iteratively checks each character. If a character in the expanded CIGAR string
#' is not an insertion ('I'), the corresponding character from the observed sequence
#' is retained in the result.
#'
#' @param cigar A CIGAR string representing the alignment of an observed sequence
#'              to a reference sequence.
#' @param obs_sequence The observed sequence (string) that is to be processed based
#'                     on the CIGAR string.
#'
#' @return A string representing the observed sequence with insertions removed.
#'
#' @keywords internal
#'
remove_insertions <- function(cigar, obs_sequence) {
  # Expanding the CIGAR string to match the length of the observed sequence
  expanded_cigar <- expand_cigar(cigar)

  # Split the observed sequence into characters
  obs_seq_split <- strsplit(obs_sequence, "")[[1]]

  # Identify the positions not marked as insertion
  non_insertion_positions <- which(substring(expanded_cigar, 1:nchar(expanded_cigar), 1:nchar(expanded_cigar)) != "I")

  # Extract the corresponding characters from the observed sequence
  result <- obs_seq_split[non_insertion_positions]

  # Combine the characters back into a single string
  return(paste(result, collapse = ""))
}
