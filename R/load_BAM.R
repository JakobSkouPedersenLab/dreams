
#' Load .BAM file into R as tibble
#'
#' @param BamPath Path to .BAM-file
#' @param chr Vector of chromosomes of interest (strings)
#' @param pos Vector of positions of interest (numeric)
#'
#' @return .BAM file in tibble format
#' @keywords internal
#'
#' @import stringr dplyr
#' @importFrom Rsamtools ScanBamParam BamFile scanBam scanBamFlag
#' @importFrom purrr map2
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#'
load_BAM <- function(BamPath, chr = NULL, pos = NULL) {



  # Get reference to BamFile
  bamFile <- BamFile(BamPath)

  # Param for loading the selected regions of BAM file
  param <- ScanBamParam(
    flag = scanBamFlag(isPaired = T, isProperPair = T, isUnmappedQuery = F, hasUnmappedMate = F, isSecondaryAlignment = F, isSupplementaryAlignment = FALSE),
    tag = c("MD", "ce", "cd", "cE", "cD"),
    which = GRanges(chr, IRanges(start = pos, end = pos)),
    what = c("qname", "rname", "strand", "pos", "mpos", "seq", "flag", "qwidth", "isize", "cigar", "mapq", "qual"),
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
    hardclip_correct_umi_features() %>%
    remove_softclips()

  return(bam_df)
}

#' Strand correct UMI features (ce and cd)
#'
#' @param df data.frame converted from lists-of-lists (scanBam)
#'
#' @return data.frame with corrected UMI features
#' @keywords internal
#'
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

#' Hard clip correct UMI features (ce and cd)
#'
#' @param df data.frame converted from lists-of-lists (scanBam)
#'
#' @return data.frame with corrected UMI features
#' @keywords internal
#'
hardclip_correct_umi_features <- function(df) {

  # If UMI features not present -> return input
  if (!all(c("ce", "cd") %in% colnames(df))) {
    return(df)
  }


  # Filter soft clips and remove from sequence string
  cigar <- df$cigar


  hard_clips_start <- ifelse(stringr::str_detect(cigar, "^[0-9]*H"),
    stringr::str_extract(string = cigar, pattern = "^[0-9]*") %>% as.numeric(),
    0
  )

  hard_clips_end <- ifelse(stringr::str_detect(cigar, "[0-9]*H$"),
    stringr::str_extract(string = cigar, pattern = "[0-9]*(?=H$)") %>% as.numeric(),
    0
  )
  trim_list <- function(x, start, end) {
    return(x[(start + 1):(length(x) - end)])
  }
  # Extract old ce and cd
  cd <- df$cd
  ce <- df$ce

  # Trim UMI features
  new_cd <- purrr::pmap(list(cd, hard_clips_start, hard_clips_end), trim_list)
  new_ce <- purrr::pmap(list(ce, hard_clips_start, hard_clips_end), trim_list)

  # Update cd and ce
  df$cd <- new_cd
  df$ce <- new_ce


  # Return updated df
  return(df)
}
