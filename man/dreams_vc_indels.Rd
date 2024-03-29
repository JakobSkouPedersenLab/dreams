% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/call_mutations_indels.R
\name{dreams_vc_indels}
\alias{dreams_vc_indels}
\title{Variant Calling in Indels}
\usage{
dreams_vc_indels(
  mutations_df,
  bam_file_path,
  reference_path,
  model,
  model_indels,
  alpha = 0.05,
  use_turboem = TRUE,
  calculate_confidence_intervals = FALSE,
  batch_size = NULL
)
}
\arguments{
\item{mutations_df}{A dataframe containing the list of mutations to be analyzed.}

\item{bam_file_path}{Path to the BAM file containing sequencing data.}

\item{reference_path}{Path to the reference genome file, typically in FASTA format.}

\item{model}{The model to be used for calling mutations for SNVs.}

\item{model_indels}{The model to be used for calling mutations for indels.}

\item{alpha}{Significance level for statistical testing, default is 0.05.}

\item{use_turboem}{Logical flag indicating whether to use the turboEM algorithm, default is TRUE.}

\item{calculate_confidence_intervals}{Logical flag indicating whether to calculate confidence intervals, default is FALSE.}

\item{batch_size}{Number of positions to process in each batch; if NULL, it's determined based on the data.}
}
\value{
A \code{\link[=data.frame]{data.frame()}} with information about the individual mutation calls, including:
\describe{
\item{chr, genomic_pos}{The genomic position of the mutation.}
\item{ref, alt}{The reference and alternative allele.}
\item{EM_converged}{If the EM algorithm converged.}
\item{EM_steps, fpeval, objfeval}{Number of steps and function evaluations by the EM algorithm.}
\item{tf_est}{The estiamted tumor fraction (allele fraction).}
\item{tf_min, tf_max}{The confidence interval of \code{tf_est}.}
\item{exp_count}{The expected count of the alternative allele under the error (null) model.}
\item{count}{The count of the alternative allele.}
\item{coverage}{The coverage used by the model (only referenceredas with and alternative allele).}
\item{full_coverage}{The total coverage of the position (for reference).}
\item{obs_freq}{The observed frequency of the alternative allele.}
\item{ll_0, ll_A}{The value of the log-likelihood function under the null (tf=0) and alternative (tf>0) hypothesis.}
\item{Q_val, df, p_val}{The chisq test statistic, degrees of freedom and p-value of the statistical test.}
\item{mutation_detected}{Whether the mutation was detected at the supplied alpha level.}
}
}
\description{
This function performs variant calling on indels from a given BAM file using a predefined list of mutations.
It processes the mutations data, batch processes the genomic positions, and calls mutations using a specific model.
}
\seealso{
\code{\link[=call_mutations]{call_mutations()}}, \code{\link[=call_cancer]{call_cancer()}}, \code{\link[=train_dreams_model]{train_dreams_model()}}
}
