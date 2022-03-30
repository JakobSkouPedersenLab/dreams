#' Get read positions from BAM-file
#'
#' @description This functions extracts information about read position that cover positions of interest in a BAM-file.
#'
#' @param bam_file_path Path to BAM-file.
#' @param chr,genomic_pos Vectors. Should specify the positions of interest (\code{chr} = Chromosome, \code{genomic_pos} = Position in chromosome)
#' @param reference_path Path to reference genome e.g. FASTA-file.
#' @param pos_wise Handle bam files position wise
#' @param chr_wise handle bam files chromosome wise
#'
#' @return [data.frame()]. Each line describes a position in a read.
#' @export
get_read_positions_from_BAM <- function(bam_file_path, chr, genomic_pos, reference_path, chr_wise = F, pos_wise = F) {

  chr_vec = chr
  pos_vec = genomic_pos

  if (chr_wise) {
    queue = lapply(unique(chr_vec), function(c) list(chr = chr_vec[chr_vec == c], pos = pos_vec[chr_vec == c]))
  } else if (pos_wise) {
    queue = lapply(1:length(chr_vec), function(i) list(chr = chr_vec[i], pos = pos_vec[i]))
  } else {
    queue = list(list(chr = chr_vec, pos = pos_vec))
  }

  read_positions_df = NULL

  for (q in queue) {

    # Load reads from BAM into data.frame
    bam_df <- load_BAM(bam_file_path, q$chr, q$pos)

    # If no coverage -> return empty data.frame
    if (nrow(bam_df) == 0) {
      return(data.frame())
    }

    # Extract features from BAM
    read_positions_df_current <-
      extract_features_from_bam(
        bam_df = bam_df,
        reference_path = reference_path
      )

    read_positions_df = rbind(read_positions_df, read_positions_df_current)

  }

  return(read_positions_df)
}
