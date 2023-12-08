#' Get read positions from BAM-file
#'
#' @description This functions extracts information about read position that cover positions of interest in a BAM-file.
#'
#' @param bam_file_path Path to BAM-file.
#' @param chr,genomic_pos Vectors. Should specify the positions of interest (\code{chr} = Chromosome, \code{genomic_pos} = Position in chromosome)
#' @param reference_path Path to reference genome e.g. FASTA-file.
#' @param batch_size Number of positions to process at a time

#'
#' @return [data.frame()]. Each line describes a position in a read.
#'
#' @keywords internal
#'
get_read_positions_from_BAM_indels <- function(bam_file_path, chr, genomic_pos, reference_path, batch_size = NULL) {

  # Only extract reads from distinct positions
  positions <- data.frame(
    chr = chr,
    genomic_pos = genomic_pos
  ) %>% distinct()

  if (is.null(batch_size)) {
    batch_size <- nrow(positions) + 1
  }

  position_batches <- positions %>% mutate(batch_idx = (row_number() %/% batch_size))

  n_batches <- length(unique(position_batches$batch_idx))


  read_positions_df <- NULL

  for (batch in sort(unique(position_batches$batch_idx))) {
    q <- position_batches %>% filter(batch_idx == batch)


    # Load reads from BAM into data.frame

    bam_df <- load_BAM(bam_file_path, q$chr, q$genomic_pos) %>%
      # Remove insertions from observed sequence based on CIGAR String
      mutate(seq_w_corr_gen_pos = mapply(remove_insertions, .data$cigar, .data$seq),
             cleaned_cigar = lapply(lapply(.data$cigar, expand_cigar), clean_insertions),
             pos_idx = .data$genomic_pos - .data$pos + 1,
             context11 =
               get_reference_seq(
                 chr = .data$chr,
                 genomic_pos = .data$genomic_pos,
                 buffer = 5,
                 reference_path = reference_path),
             ref = substring(.data$context11, 6, 6),) %>%
               mutate(obs = case_when(
                 substring(cleaned_cigar, pos_idx, pos_idx) == "D" ~ "D",
                 substring(cleaned_cigar, pos_idx, pos_idx) == "I" ~ "I",
                 substring(cleaned_cigar, pos_idx, pos_idx) == "M" &
                 substring(seq, pos_idx, pos_idx) == ref ~ substring(seq, pos_idx, pos_idx),
                 TRUE ~ "N"
               ))



    # If no coverage -> return empty data.frame
    if (nrow(bam_df) == 0) {
      return(data.frame())
    }

    # Extract features from BAM
    read_positions_df_current <- bam_df

    read_positions_df <- rbind(read_positions_df, read_positions_df_current)
  }

  return(read_positions_df)
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
