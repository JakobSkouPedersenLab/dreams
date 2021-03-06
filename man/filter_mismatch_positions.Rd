% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_training_data.R
\name{filter_mismatch_positions}
\alias{filter_mismatch_positions}
\title{Title}
\usage{
filter_mismatch_positions(
  read_positions,
  bam_file,
  mm_rate_max = 1,
  bed_include_path = NULL,
  positions_to_exclude_paths = NULL
)
}
\arguments{
\item{read_positions}{dataframe of read positions}

\item{bam_file}{bam file path}

\item{mm_rate_max}{maximum mm_rate for positions}

\item{bed_include_path}{bed regions to include in training data}

\item{positions_to_exclude_paths}{positions to exclude from training}
}
\value{
filtered read position dataframe
}
\description{
Title
}
\keyword{internal}
