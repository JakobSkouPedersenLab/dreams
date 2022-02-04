#' Load .BAM file into R as tibble
#'
#' @param BamPath Path to .BAM-file
#' @param chr Vector of chromosomes of interest (strings)
#' @param pos Vector of positions of interest (numeric)
#'
#' @return .BAM file in tibble format
#'
#' @import Rsamtools stringr dplyr
#' @importFrom purrr map2
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
load_BAM <- function(BamPath, chr = NULL, pos = NULL) {

  # Function to read bam file, and perform initial filtering

  bamFile <- BamFile(BamPath)
  print(paste0("PATH: ", BamPath))

  if (!is.null(chr)) {
    # Param for loading the selected regions of BAM file
    param <- ScanBamParam(
      tag = c("MD", "ce", "cd", "cE", "cD"),
      which = GRanges(chr, IRanges(start = pos, end = pos)),
      what = c("qname", "rname", "strand", "pos", "mpos", "seq", "flag", "qwidth", "isize", "cigar", "mapq", "qual")
    )
  } else {
    # Param for loading the entire BAM file
    param <- ScanBamParam(
      tag = c("MD", "ce", "cd", "cE", "cD"),
      what = c("qname", "rname", "strand", "pos", "mpos", "seq", "flag", "qwidth", "isize", "cigar", "mapq", "qual")
    )
  }

  bam <- scanBam(bamFile, param = param)

  for (i in 1:length(bam)) {
    # Unpack MD tag
    bam[[i]]$MD <- str_to_upper(bam[[i]]$tag$MD)
    bam[[i]]$ce <- bam[[i]]$tag$ce
    bam[[i]]$cd <- bam[[i]]$tag$cd
    bam[[i]]$cE <- bam[[i]]$tag$cE
    bam[[i]]$cD <- bam[[i]]$tag$cD
    bam[[i]]$tag <- NULL
    bam[[i]]$seq <- as.character(bam[[i]]$seq)
    bam[[i]]$qual <- as.character(bam[[i]]$qual)
    bam[[i]]$chr <- as.character(bam[[i]]$rname)
    # Make genomic position into features

    if (!is.null(chr)) {
      genomic_pos <- as.numeric(as.character(str_extract(names(bam[i]), "[0-9]*$")))
      genomic_pos_vec <- rep(genomic_pos, length(bam[[i]]$rname))
      bam[[i]]$genomic_pos <- genomic_pos_vec
    }
  }

  # TODO: Look at optimization
  # Possible speed optimization:
  # https://stackoverflow.com/questions/44652234/read-bam-with-tags-to-tbl-recursive-dplyrbind-cols-for-list-of-lists

  # Convert list-of-lists into data.frame + make features and names nicer

  bam_df <-
    bind_rows(lapply(bam, as_tibble)) %>%
    #    rename(
    #     chr = rname
    #  ) %>%
    mutate(
      strand =
        ifelse(
          as.character(.data$strand) == "+",
          "fwd",
          "rev"
        ),
      chr = as.character(chr)
    )

  print(dim(bam_df))


  # TODO: Move filtering to other function!
  # Initial filtering
  if (nrow(bam_df) == 0) {
    return(bam_df)
  }

  # print(head(bam_df))

  bam_df <- bam_df %>%
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
    ) %>%
    remove_softclips()

  return(bam_df)
}
