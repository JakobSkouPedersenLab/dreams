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
#' @export
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

    bam_df <- load_BAM(bam_file_path, q$chr, q$genomic_pos)


    # If no coverage -> return empty data.frame
    if (nrow(bam_df) == 0) {
      return(data.frame())
    }

    # Extract features from BAM
    read_positions_df_current <-
      extract_features_from_bam_indels(
        bam_df = bam_df,
        reference_path = reference_path
      )

    read_positions_df <- rbind(read_positions_df, read_positions_df_current)
  }

  return(read_positions_df)
}
