
#' Load .BAM file into R as tibble
#'
#' @param BamPath Path to .BAM-file
#' @param chr Vector of chromosomes of interest (strings)
#' @param pos Vector of positions of interest (numeric)
#'
#' @return .BAM file in tibble format
#'
#' @import stringr dplyr
#' @importFrom Rsamtools ScanBamParam BamFile scanBam
#' @importFrom purrr map2
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
load_BAM <- function(BamPath, chr = NULL, pos = NULL) {

  # Get reference to BamFile
  bamFile <- BamFile(BamPath)

  # Param for loading the selected regions of BAM file
  param <- ScanBamParam(
    tag = c("MD", "ce", "cd", "cE", "cD"),
    which = GRanges(chr, IRanges(start = pos, end = pos)),
    what = c("qname", "rname", "strand", "pos", "mpos", "seq", "flag", "qwidth", "isize", "cigar", "mapq", "qual")
  )

  # Load BAM file
  bam <- scanBam(bamFile, param = param)

  # Unpack tags
  for (i in 1:length(bam)) {
    # Default tags
    bam[[i]]$MD <- str_to_upper(bam[[i]]$tag$MD)
    bam[[i]]$seq <- as.character(bam[[i]]$seq)
    bam[[i]]$qual <- as.character(bam[[i]]$qual)
    bam[[i]]$chr <- as.character(bam[[i]]$rname)

    # UMI tags
    umi_is_present <- all(!sapply(bam[[i]]$tag[c("ce", "cd", "cE", "cD")], is.null))
    if (umi_is_present) {
      bam[[i]]$ce <- bam[[i]]$tag$ce
      bam[[i]]$cd <- bam[[i]]$tag$cd
      bam[[i]]$cE <- bam[[i]]$tag$cE
      bam[[i]]$cD <- bam[[i]]$tag$cD
    }

    bam[[i]]$tag <- NULL

    # Make genomic position into feature
    if (!is.null(chr) & !is.null(pos)) {
      genomic_pos <- as.numeric(as.character(str_extract(names(bam[i]), "[0-9]*$")))
      bam[[i]]$genomic_pos <- rep(genomic_pos, length(bam[[i]]$rname))
    }
  }

  # Convert list-of-lists into data.frame + make features and names nicer
  raw_bam_df <- bind_rows(lapply(bam, as_tibble))

  # If no reads -> Return data.frame
  if (nrow(raw_bam_df) == 0) {
    return(raw_bam_df)
  }

  # Clean up columns + Reverse ce and cd tag on reverse strand
  bam_df <- raw_bam_df %>%
    mutate(
      strand = ifelse(
        as.character(.data$strand) == "+", "fwd", "rev"
      )
    ) %>%
    strand_correct_umi_features() %>%
    remove_softclips()

  return(bam_df)
}

#' Strand correct UMI features (ce and cd)
#'
#' @param df data.frame converted from lists-of-lists (scanBam)
#'
#' @return data.frame with corrected UMI features
strand_correct_umi_features <- function(df) {

  # If UMI features not present -> return input
  if (!all(c("ce", "cd") %in% colnames(df))) {
    return(df)
  }

  df %>%
    mutate(
      ce = map2(
        .data$strand,
        .data$ce,
        function(strand, ce) {
          if (strand == "rev") rev(ce) else ce
        }
      ),
      cd = map2(
        .data$strand,
        .data$cd,
        function(strand, cd) {
          if (strand == "rev") rev(cd) else cd
        }
      )
    )
}
