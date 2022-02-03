

#' Remove softclips
#'
#' @param df A dataframe from load_BAM
#'
#' @return A new dataframe with softclipped bases removed
#'
remove_softclips <- function(df) {
  # Extract seq and cigar
  seq <- df$seq
  qual <- df$qual
  cigar <- df$cigar
  cd <- df$cd
  ce <- df$ce

  # Filter soft clips and remove from sequence string
  soft_clips_start <- ifelse(stringr::str_detect(cigar, "^[0-9]*S"),
    stringr::str_extract(string = cigar, pattern = "^[0-9]*") %>% as.numeric(),
    0
  )

  soft_clips_end <- ifelse(stringr::str_detect(cigar, "[0-9]*S$"),
    stringr::str_extract(string = cigar, pattern = "[0-9]*(?=S$)") %>% as.numeric(),
    0
  )

  trim_list <- function(x, start, end) {
    return(x[(start + 1):(length(x) - end)])
  }

  # Trim cd
  new_cd <- purrr::pmap(list(cd, soft_clips_start, soft_clips_end), trim_list)

  # Trim ce
  new_ce <- purrr::pmap(list(ce, soft_clips_start, soft_clips_end), trim_list)

  # Trim sequence
  new_seq <- substring(seq, soft_clips_start + 1, nchar(seq) - soft_clips_end)

  # Trim qual
  new_qual <- substring(qual, soft_clips_start + 1, nchar(qual) - soft_clips_end)

  # Remove soft clipped ends from cigar
  new_cigar <- stringr::str_remove_all(cigar, "^[0-9]*S|[0-9]*S$")

  # Update seq and cigar
  df$seq <- new_seq
  df$cigar <- new_cigar
  df$qual <- new_qual
  df$cd <- new_cd
  df$ce <- new_ce

  return(df)
}
