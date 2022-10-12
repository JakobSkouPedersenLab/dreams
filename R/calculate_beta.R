#' Calculate beta-factor
#'
#' @description This functions calculates the beta factor based on a BAM-file.
#' @param bam_file_path Path to .BAM-file
#' @param factor mismatch to match ratio in training data
#' @param bed_file path to bed file with relevant regions
#' @param mm_rate_max maximum error rate for beta calculation
#' @param reference_path Path to reference genome e.g. FASTA-file.
#'
#' @return [data.frame()]. Each line describes a position in a read.
#'
#' @export
calculate_beta_factor <- function(bam_file_path, factor, reference_path, bed_file = NULL, mm_rate_max = NULL) {
  pp <- Rsamtools::PileupParam(
    max_depth = 250000000, min_base_quality = 13, min_mapq = 0,
    min_nucleotide_depth = 1, min_minor_allele_depth = 0,
    distinguish_strands = FALSE, distinguish_nucleotides = T,
    ignore_query_Ns = TRUE, include_deletions = TRUE, include_insertions = FALSE,
    left_bins = NULL, query_bins = NULL, cycle_bins = NULL
  )



  if (!is.null(bed_file)) {
    included_regions_granges <- bed_to_granges(bed_file)
    coverage_data <- Rsamtools::pileup(bam_file_path, pileupParam = pp, scanBamParam = ScanBamParam(which = included_regions_granges))
  } else {
    coverage_data <- Rsamtools::pileup(bam_file_path, pileupParam = pp)
  }

  coverage_data <- coverage_data %>%
    rename(
      chr = seqnames,
      genomic_pos = pos,
      coverage = count
    )

  coverage_data <- coverage_data %>%
    mutate(reference_base = get_reference_seq(.data$chr, .data$genomic_pos, 0, reference_path))

  coverage_data <- coverage_data %>%
    group_by(.data$chr, .data$genomic_pos) %>%
    mutate(total_coverage = sum(.data$coverage)) %>%
    ungroup()

  errors <- coverage_data %>% filter(
    .data$reference_base != .data$nucleotide,
    .data$coverage < mm_rate_max * .data$total_coverage
  )


  beta <- sum(errors$coverage) / (sum(coverage_data$coverage) - sum(errors$coverage))


  return(beta)
}
