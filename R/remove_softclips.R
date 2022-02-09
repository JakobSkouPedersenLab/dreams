

#' Remove softclips
#'
#' @param df A dataframe from load_BAM
#'
#' @return A new dataframe with softclipped bases removed
#'
remove_softclips <- function(df) {
  # Extract old seq, qual and cigar
  seq <- df$seq
  qual <- df$qual
  cigar <- df$cigar

  # Filter soft clips and remove from sequence string
  soft_clips_start <- ifelse(stringr::str_detect(cigar, "^[0-9]*S"),
    stringr::str_extract(string = cigar, pattern = "^[0-9]*") %>% as.numeric(),
    0
  )

  soft_clips_end <- ifelse(stringr::str_detect(cigar, "[0-9]*S$"),
    stringr::str_extract(string = cigar, pattern = "[0-9]*(?=S$)") %>% as.numeric(),
    0
  )

  # Trim sequence, qual and cigar
  new_seq <- substring(seq, soft_clips_start + 1, nchar(seq) - soft_clips_end)
  new_qual <- substring(qual, soft_clips_start + 1, nchar(qual) - soft_clips_end)
  new_cigar <- stringr::str_remove_all(cigar, "^[0-9]*S|[0-9]*S$")

  # Update seq and cigar
  df$seq <- new_seq
  df$cigar <- new_cigar
  df$qual <- new_qual

  # Handle UMI features if present
  if (all(c("ce", "cd") %in% colnames(df))) {
    trim_list <- function(x, start, end) {
      return(x[(start + 1):(length(x) - end)])
    }
    # Extract old ce and cd
    cd <- df$cd
    ce <- df$ce

    # Trim UMI features
    new_cd <- purrr::pmap(list(cd, soft_clips_start, soft_clips_end), trim_list)
    new_ce <- purrr::pmap(list(ce, soft_clips_start, soft_clips_end), trim_list)

    # Update cd and ce
    df$cd <- new_cd
    df$ce <- new_ce
  }

  # Return updated df
  return(df)
}
