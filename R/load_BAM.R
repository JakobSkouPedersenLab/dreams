
#' Load .BAM file into R as tibble
#'
#' @description This function reads alignment data from a .BAM file and formats it as a tibble (data frame) in R.
#' It allows for the selection of specific genomic regions by specifying chromosomes and positions, and it extracts
#' relevant alignment data including tags for further analysis. The function also includes processing of UMI tags
#' and correction for hard clips in the CIGAR strings.
#'
#' @param BamPath The file path to the .BAM file to be loaded.
#' @param chr A vector of chromosome identifiers indicating the chromosomes of interest for the analysis. (strings)
#' @param pos A vector of numeric positions indicating the specific locations of interest within the chromosomes. (numeric)
#'
#' @return A tibble (data frame) containing the extracted alignment data from the .BAM file, including processed UMI
#' features and genomic positions, as well as any necessary corrections for alignment issues such as hard clipping.
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

  # Set up the parameters for scanning the BAM file, focusing on paired and properly paired alignments,
  # excluding unmapped queries, secondary alignments, and supplementary alignments.
  # Tags 'MD', 'ce', 'cd', 'cE', and 'cD' are included for extraction.
  param <- ScanBamParam(
    flag = scanBamFlag(isPaired = T, isProperPair = T, isUnmappedQuery = F, hasUnmappedMate = F, isSecondaryAlignment = F, isSupplementaryAlignment = FALSE),
    tag = c("MD", "ce", "cd", "cE", "cD"),
    which = GRanges(chr, IRanges(start = pos, end = pos)),
    what = c("qname", "rname", "strand", "pos", "mpos", "seq", "flag", "qwidth", "isize", "cigar", "mapq", "qual"),
  )

  # Load the BAM file based on the parameters set above using the scanBam function from Rsamtools.
  bam <- scanBam(bamFile, param = param)


  # Process each alignment entry, unpacking and formatting the extracted data.
  for (i in 1:length(bam)) {
    # Convert MD tag to upper case, sequence and quality scores to character vectors, and the reference name as well.
    bam[[i]]$MD <- str_to_upper(bam[[i]]$tag$MD)
    bam[[i]]$seq <- as.character(bam[[i]]$seq)
    bam[[i]]$qual <- as.character(bam[[i]]$qual)
    bam[[i]]$chr <- as.character(bam[[i]]$rname)

    # Check for presence of UMI tags and if present, add them to the list.
    umi_is_present <- all(!sapply(bam[[i]]$tag[c("ce", "cd", "cE", "cD")], is.null))
    if (umi_is_present) {
      bam[[i]]$ce <- bam[[i]]$tag$ce
      bam[[i]]$cd <- bam[[i]]$tag$cd
      bam[[i]]$cE <- bam[[i]]$tag$cE
      bam[[i]]$cD <- bam[[i]]$tag$cD
    }
    # Remove the tag list to tidy up the data structure.
    bam[[i]]$tag <- NULL

    # If chromosomes and positions are specified, extract genomic positions from the alignment names.
    if (!is.null(chr) & !is.null(pos)) {
      genomic_pos <- as.numeric(as.character(str_extract(names(bam[i]), "[0-9]*$")))
      bam[[i]]$genomic_pos <- rep(genomic_pos, length(bam[[i]]$rname))
    }

  }

  # Combine the list of tibbles into one data frame, tidying up the column names in the process.
  raw_bam_df <- bind_rows(lapply(bam, as_tibble))

  # If the dataframe is empty, return it as is.
  if (nrow(raw_bam_df) == 0) {
    return(raw_bam_df)
  }

  # Further process the dataframe by correcting the strand information for UMI features,
  # correcting UMI features for hard clips, and removing soft clips from the data.
  bam_df <- raw_bam_df %>%
    mutate(
      strand = ifelse(
        as.character(.data$strand) == "+", "fwd", "rev"
      )
    ) %>%
    strand_correct_umi_features() %>%
    hardclip_correct_umi_features() %>%
    remove_softclips() %>%
    remove_hardclips()

  # Return the fully processed tibble.
  return(bam_df)
}

#' Strand correct UMI features (ce and cd)
#'
#' @description This function corrects the orientation of Unique Molecular Identifier (UMI) features based on
#' the DNA strand from which they originate. It takes a data frame containing UMI features and adjusts them
#' if they are from the reverse strand, ensuring that the orientation of the UMI sequences is consistent with
#' their representation on the forward strand. This correction is necessary for accurate downstream analysis
#' of sequence data.
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

  # Use dplyr's mutate to apply changes to the ce and cd columns
  df %>%
    mutate(
      # Use purrr's map2 function to iterate over elements of .data$strand and .data$ce together
      ce = map2(
        .data$strand,
        .data$ce,
        function(strand, ce) {
          # If the strand is reverse, reverse the elements of the ce array
          if (strand == "rev") rev(ce) else ce
        }
      ),
      # Do the same for the cd column
      cd = map2(
        .data$strand,
        .data$cd,
        function(strand, cd) {
          # If the strand is reverse, reverse the elements of the cd array
          if (strand == "rev") rev(cd) else cd
        }
      )
    )
}

#' Hard clip correct UMI features (ce and cd)
#'
#' @description This function identifies and corrects for hard clipped bases in the UMI (Unique Molecular Identifier)
#' features within sequencing data. Hard clipping, indicated by 'H' in the CIGAR string of a BAM file, refers to
#' bases of the sequence that are not aligned to the reference genome and are therefore not included in the read.
#' This function adjusts the UMI feature sequences accordingly to ensure that only the aligned portions of the UMIs
#' are considered in downstream analyses. The function applies corrections to both 'ce' and 'cd' UMI feature columns
#' in the provided dataframe.
#'
#' @param df data.frame converted from lists-of-lists (scanBam)
#'
#' @return data.frame with corrected UMI features
#' @keywords internal
#'
hardclip_correct_umi_features <- function(df) {

  # Check if the UMI feature columns ('ce' and 'cd') are present in the dataframe; if not, return the dataframe unmodified
  if (!all(c("ce", "cd") %in% colnames(df))) {
    return(df)
  }


  # Extract the CIGAR string for each read from the dataframe
  cigar <- df$cigar

  # Calculate the number of bases hard clipped at the start of each read
  hard_clips_start <- ifelse(stringr::str_detect(cigar, "^[0-9]*H"),
    stringr::str_extract(string = cigar, pattern = "^[0-9]*") %>% as.numeric(),
    0
  )
  # Calculate the number of bases hard clipped at the end of each read
  hard_clips_end <- ifelse(stringr::str_detect(cigar, "[0-9]*H$"),
    stringr::str_extract(string = cigar, pattern = "[0-9]*(?=H$)") %>% as.numeric(),
    0
  )
  # Define a helper function to trim a vector based on the start and end positions
  trim_list <- function(x, start, end) {
    return(x[(start + 1):(length(x) - end)])
  }
  # Extract the original UMI features 'cd' and 'ce'
  cd <- df$cd
  ce <- df$ce

  # Use purrr's pmap function to apply the trimming to each UMI feature
  new_cd <- purrr::pmap(list(cd, hard_clips_start, hard_clips_end), trim_list)
  new_ce <- purrr::pmap(list(ce, hard_clips_start, hard_clips_end), trim_list)

  # Update the 'cd' and 'ce' columns in the dataframe with the trimmed UMI features
  df$cd <- new_cd
  df$ce <- new_ce


  # Return the updated dataframe
  return(df)
}


#' Remove softclips
#'
#' @description This function processes a dataframe obtained from a BAM file and removes soft clipped bases from the sequences.
#' Soft clipping (indicated by 'S' in the CIGAR string) refers to bases that are present in the sequencing read but not used
#' in the alignment to the reference genome. This function trims these soft clipped bases from the sequence (`seq`) and quality
#' scores (`qual`), and also adjusts the CIGAR string accordingly. If UMI features are present, it will also trim these.
#'
#' @param df A dataframe from load_BAM
#' @keywords internal

#' @return A new dataframe with softclipped bases removed
#'
remove_softclips <- function(df) {
  # Extract sequence, quality scores, and CIGAR strings from the dataframe
  seq <- df$seq
  qual <- df$qual
  cigar <- df$cigar

  # Identify and calculate the length of soft clipping at the start of the sequences
  soft_clips_start <- ifelse(stringr::str_detect(cigar, "^[0-9]*S"),
                             stringr::str_extract(string = cigar, pattern = "^[0-9]*") %>% as.numeric(),
                             0
  )

  soft_clips_end <- ifelse(stringr::str_detect(cigar, "[0-9]*S$"),
                           stringr::str_extract(string = cigar, pattern = "[0-9]*(?=S$)") %>% as.numeric(),
                           0
  )

  # Trim the sequence and quality strings to remove soft clipped bases
  new_seq <- substring(seq, soft_clips_start + 1, nchar(seq) - soft_clips_end)
  new_qual <- substring(qual, soft_clips_start + 1, nchar(qual) - soft_clips_end)

  # Remove soft clipping indications from the CIGAR strings
  new_cigar <- stringr::str_remove_all(cigar, "^[0-9]*S|[0-9]*S$")

  # Update the dataframe with the trimmed sequence and CIGAR strings
  df$seq <- new_seq
  df$cigar <- new_cigar
  df$qual <- new_qual

  # If UMI feature columns are present, apply trimming to these as well
  if (all(c("ce", "cd") %in% colnames(df))) {
    trim_list <- function(x, start, end) {
      return(x[(start + 1):(length(x) - end)])
    }
    # Extract the current UMI feature values
    cd <- df$cd
    ce <- df$ce

    # Trim UMI features for soft clipped bases
    new_cd <- purrr::pmap(list(cd, soft_clips_start, soft_clips_end), trim_list)
    new_ce <- purrr::pmap(list(ce, soft_clips_start, soft_clips_end), trim_list)

    # Update the dataframe with the trimmed UMI features
    df$cd <- new_cd
    df$ce <- new_ce
  }

  # Return the dataframe with soft clipped bases removed
  return(df)
}


#' Remove Hard Clipping from CIGAR Strings in a DataFrame
#'
#' @description This function processes a dataframe containing CIGAR strings
#' in one of its columns. It removes the hard clipping ('H') and deletions ('D')
#' at the beginning and end of each CIGAR string. Hard clipping and deletions
#' do not contribute to the alignment sequence, and their removal simplifies
#' further processing of these strings.
#'
#' @param df A dataframe with a column named 'cigar', which contains CIGAR strings.
#'
#' @return A dataframe identical to the input, except the CIGAR strings in the
#' 'cigar' column are modified to exclude initial and terminal hard clipping
#' and deletions.
#'
#'@keywords internal
remove_hardclips <- function(df){
  cigar <- df$cigar
  cigar <- str_remove_all(cigar, pattern = "^[0-9]+[HD]")
  cigar <- str_remove_all(cigar, pattern = "[0-9]+[HD]$")
  df$cigar = cigar
  return(df)
}



