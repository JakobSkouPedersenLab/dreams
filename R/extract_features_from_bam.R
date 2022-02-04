#' Title
#'
#' @param s
#' @param k
#'
#' @return
#'
#' @importFrom stringi stri_count_fixed
calc_string_entropy_k_mer <- function(s, k = 2) {

  # Generate k-mers
  alphabet <- c("A", "C", "G", "T", "N")
  alphabet_k_rep_list <- rep(list(alphabet), k)
  k_mer_df <- expand.grid(alphabet_k_rep_list)
  k_mer_vec <- apply(k_mer_df, 1, paste0, collapse = "")

  s_length <- nchar(s) - (k - 1)

  f_mat <- t(sapply(s, stri_count_fixed, pattern = k_mer_vec, overlap = TRUE)) / s_length

  # Shannon entropy
  H <- -rowSums(f_mat * log10(f_mat), na.rm = TRUE)

  return(H)
}

#' Title
#'
#' @param bam_df dataframe from load_BAM
#' @param reference_path reference genome path
#'
#' @return dataframe with read positions
#'
#' @importFrom purrr map2_int
extract_features_from_bam <- function(bam_df, reference_path) {

  # Make genomic position features
  genomic_pos_feature_df <-
    distinct(
      data.frame(
        chr = bam_df$chr,
        genomic_pos = bam_df$genomic_pos,
        strand = bam_df$strand
      )
    ) %>%
    mutate(
      context11 = get_reference_seq(
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
      trinucleotide_ctx = paste0(.data$ctx_minus1, .data$ref, .data$ctx_plus1),
      trinucleotide_ctx_strand = paste0(.data$ctx_minus1, .data$ref, .data$ctx_plus1, "_", .data$strand)
    )

  # Make read specific features
  read_feature_df <-
    bam_df %>%
    mutate(
      pos_idx = .data$genomic_pos - .data$pos + 1
    ) %>%
    correct_pos_idx_w_cigar() %>%
    mutate(
      obs = ifelse(
        .data$is_in_deletion,
        "N",
        substring(.data$seq, .data$pos_idx, .data$pos_idx)
      ),
      fragment_size = abs(.data$isize),
      seq_length = nchar(.data$seq),
      read_index = if_else(.data$strand == "fwd", .data$pos_idx, .data$seq_length - .data$pos_idx + 1),
      first_in_pair = as.integer(as.logical(bitwAnd(.data$flag, 64))),
      umi_count = map2_int(.data$cd, .data$pos_idx, function(cd, pos_idx) cd[pos_idx]),
      umi_errors = map2_int(.data$ce, .data$pos_idx, function(ce, pos_idx) ce[pos_idx]),
      n_errors_in_read = str_count(.data$MD, "\\d+[ATCG]"),
      n_insertions_in_read = str_count(.data$cigar, "I"),
      n_deletions_in_read = str_count(.data$cigar, "D")
    ) %>%
    # TODO: Move to filter function! Or do before calling this function!
    filter(.data$fragment_size != 0)

  # Join and select features: Read, genomic positions and UMI
  feature_df <-
    left_join(read_feature_df, genomic_pos_feature_df,
      by = c("chr", "genomic_pos", "strand")
    ) %>%
    mutate(n_other_errors = .data$n_errors_in_read - ifelse((!.data$is_in_deletion) & (.data$obs != .data$ref), 1, 0))

  feature_df <- feature_df %>%
    select(
      .data$qname, .data$chr, .data$genomic_pos, .data$obs, .data$ref,
      strand, .data$first_in_pair, .data$read_index, .data$fragment_size,
      ctx_minus1, .data$ctx_plus1, .data$trinucleotide_ctx, .data$trinucleotide_ctx_strand, .data$context11,
      local_complexity_1, .data$local_complexity_2, .data$local_GC, .data$umi_count, .data$umi_errors,
      n_other_errors, .data$n_insertions_in_read, .data$n_deletions_in_read
    )
  return(feature_df)
}

#' Get reference sequence
#'
#' @param chr chromosome
#' @param genomic_pos genomic_position
#' @param buffer how many neighbors to include
#' @param reference_path reference genome fa
#'
#' @importFrom Rsamtools FaFile getSeq

get_reference_seq <- function(chr, genomic_pos, buffer, reference_path) {
  FaFile <- FaFile(reference_path)

  Fa <- getSeq(FaFile, param = GRanges(chr, IRanges(genomic_pos - buffer, genomic_pos + buffer)))

  return(as.character(Fa))
}
