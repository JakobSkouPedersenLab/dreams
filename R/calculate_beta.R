#' calculate beta
#'
#' @export

predict_error_rates <- function(bam_file_path, factor, bed_file = NULL, mm_rate_max = NULL) {
  pp <- Rsamtools::PileupParam(
    max_depth = 250000000, min_base_quality = 13, min_mapq = 0,
    min_nucleotide_depth = 1, min_minor_allele_depth = 0,
    distinguish_strands = FALSE, distinguish_nucleotides = T,
    ignore_query_Ns = TRUE, include_deletions = TRUE, include_insertions = FALSE,
    left_bins = NULL, query_bins = NULL, cycle_bins = NULL
  )


  if (!is.null(bed_file)) {
    included_regions_granges <- bed_to_granges(bed_include_path)
    coverage_data <- Rsamtools::pileup(bam_file, pileupParam = pp, scanBamParam = ScanBamParam(which = included_regions_granges))
  } else {
    coverage_data <- Rsamtools::pileup(bam_file, pileupParam = pp)
  }

  coverage_data <- coverage_data %>%
    dplyr::rename(chr = seqnames, genomic_pos = pos, coverage = count) %>%
    mutate(reference_base = get_reference_seq(chr, genomic_pos, 0, reference_path)) %>%
    group_by(chr, genomic_pos) %>%
    mutate(total_coverage = sum(coverage)) %>%
    ungroup()

  errors <- coverage_data %>% filter(reference_base != nucleotide, coverage < mm_rate_max * total_coverage)

  beta <- sum(errors$coverage) / (sum(coverage_data$coverage) - sum(errors$coverage))


  return(beta)
}
