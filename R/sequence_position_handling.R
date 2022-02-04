#' Title correct_pos_idx_w_cigar
#'
#' @param df dataframe_with_bam_data
#' @return df
#' @importFrom rlang .data
#' @importFrom purrr map pmap pmap_dbl pmap_lgl
#' @importFrom stringr str_extract_all str_remove_all str_locate_all
correct_pos_idx_w_cigar <- function(df) {
  df %>% mutate(
    # Raw processing of CIGAR
    cigar_pos = str_extract_all(string = cigar, pattern = "\\d+(?=[MID]+)") %>% lapply(as.numeric),
    insert_idx = str_remove_all(string = cigar, pattern = "[\\d]") %>%
      str_locate_all(pattern = "I") %>%
      map(function(x) x[, "start"]),
    del_idx = .data$cigar %>%
      str_remove_all(pattern = "[\\d]") %>%
      str_locate_all(pattern = "D") %>%
      map(function(x) x[, "start"]),
    ind_sz = map2(cigar_pos, insert_idx, function(x, y) x[y]),
    del_sz = map2(cigar_pos, del_idx, function(x, y) x[y]),
    insert_read_index =
      map2(
        cigar_pos, insert_idx,
        function(cigar_pos, insert_idx) {
          cigar_pos[insert_idx] <- 0
          cumsum(cigar_pos)[insert_idx]
        }
      ),
    del_read_index = pmap(
      list(
        cigar_pos,
        del_idx,
        insert_idx
      ),
      function(cigar_pos, del_idx, insert_idx) {
        cigar_pos[insert_idx] <- 0
        cumsum(cigar_pos)[del_idx]
      }
    ),
    # Processing of read position
    ind_cor = pmap_dbl(
      list(
        ind_sz,
        insert_read_index,
        pos_idx
      ),
      function(ind_sz, insert_read_index, pos_idx) {
        sum(ind_sz[insert_read_index < pos_idx])
      }
    ),
    del_cor = pmap_dbl(
      list(
        del_sz,
        del_read_index,
        pos_idx
      ),
      function(del_sz, del_read_index, pos_idx) {
        sum(del_sz[del_read_index < pos_idx])
      }
    ),
    is_in_deletion = pmap_lgl(
      list(
        del_sz,
        del_read_index,
        pos_idx
      ),
      function(del_sz, del_read_index, pos_idx) {
        any((del_read_index - del_sz + 1 <= pos_idx) & (pos_idx <= del_read_index))
      }
    ),
    old_pos_idx = pos_idx,
    pos_idx = pos_idx + ind_cor - del_cor
  )
}


#' Title get_mismatch_genomic_pos_list
#'
#' @param pos
#' @param MDtag
#'
#' @return genomic positions
#' @importFrom stringr str_extract_all
#' @importFrom stringr str_count
#' @importFrom stringr str_which
#' @importFrom purrr map2
#' @importFrom rlang .data
get_mismatch_genomic_pos_list <- function(pos, MDtag) {
  # Number of leading matches
  leading_matches <- str_extract_all(MDtag, "\\d+(?=[\\^ATCGN]+)") %>% lapply(as.numeric)

  event_length <- str_extract_all(MDtag, "[\\^ATCGN]+") %>% lapply(str_count, pattern = "[ATCGN]")

  pos_idx_list <- map2(leading_matches, event_length, function(x, y) cumsum(x + y))

  mismatch_idx <-
    str_extract_all(MDtag, "[\\^ATCGN]+") %>%
    lapply(str_which, pattern = "\\^", negate = TRUE)

  mismatch_genomic_positions_list <-
    map2(pos_idx_list, mismatch_idx, function(x, y) x[y]) %>%
    map2(pos, function(x, y) x + y - 1)
  return(mismatch_genomic_positions_list)
}




#' Title get_match_genomic_pos_list
#'
#' @param pos position_list
#' @param cigar cigar string list
#' @param MDtag mdtag list
#'
#' @return match_positions
get_match_genomic_pos_list <- function(pos, cigar, MDtag) {
  # Remove insert from cigar
  cigar_inserts_removed <- cigar %>% str_remove_all(pattern = "[\\d]+I")

  event_lengths <-
    cigar_inserts_removed %>%
    str_extract_all(pattern = "\\d+(?=[MD]+)") %>%
    lapply(as.numeric)

  genomic_offset <- map(event_lengths, function(x) 1:sum(x))

  del_idx <-
    str_extract_all(cigar_inserts_removed, "[MD]") %>%
    lapply(str_which, pattern = "D")

  del_start <-
    map(event_lengths, cumsum) %>%
    map2(event_lengths, function(x, y) x - y + 1) %>%
    map2(del_idx, function(x, y) x[y])

  del_end <-
    map(event_lengths, cumsum) %>%
    map2(del_idx, function(x, y) x[y])

  del_positions <- map2(del_start, del_end, function(x, y) mapply(":", x, y) %>% unlist())
  mismatch_genomic_pos <- get_mismatch_genomic_pos_list(pos = pos, MDtag = MDtag)

  match_positions <-
    genomic_offset %>%
    # Remove deletions
    map2(del_positions, function(x, y) setdiff(x, y)) %>%
    # Get genomic position by adding read start position
    map2(pos, function(x, y) x + y - 1) %>%
    # Filter mismatch positions
    map2(mismatch_genomic_pos, function(x, y) setdiff(x, y))

  return(match_positions)
}
