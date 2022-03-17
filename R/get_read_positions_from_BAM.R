#' Get read positions from BAM-file
#'
#' @description This functions extracts information about read position that cover positions of interest in a BAM-file.
#'
#' @param bam_file_path Path to BAM-file.
#' @param chr,genomic_pos Vectors. Should specify the positions of interest (\code{chr} = Chromosome, \code{genomic_pos} = Position in chromosome)
#' @param reference_path Path to reference genome e.g. FASTA-file.
#'
#' @return [data.frame()]. Each line describes a position in a read.
#' @export
get_read_positions_from_BAM <- function(bam_file_path, chr, genomic_pos, reference_path) {

  # Load reads from BAM into data.frame
  bam_df <- load_BAM(bam_file_path, chr, genomic_pos)

  # If no coverage -> return empty data.frame
  if (nrow(bam_df) == 0) {
    return(data.frame())
  }

  # Extract features from BAM
  read_positions_df <-
    extract_features_from_bam(
      bam_df = bam_df,
      reference_path = reference_path
    )

  return(read_positions_df)
}
