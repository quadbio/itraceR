% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_merge_barcodes_scars.r
\name{do_filter_barcodes}
\alias{do_filter_barcodes}
\title{Filter the iTracer barcode data frame}
\usage{
do_filter_barcodes(
  df_barcode,
  coverage_cutoff_mode = c("linear", "log_gmm", "log_turn"),
  cutoff_prob = 1 - 1e-10,
  num_log_breaks = 20,
  do_plot = FALSE
)
}
\arguments{
\item{df_barcode}{Barcode data frame}

\item{coverage_cutoff_mode}{Algorithm to determine read coverage threshold. 'linear' is the default method described in the paper. 'log_gmm' applies a bimodal Gaussian mixture model. 'log_turn' is similar to 'linear' but with read number log-transformed.}

\item{cutoff_prob}{Cutoff of the bimodal Gaussian mixture model posterior probability. Only applicable when coverage_cutoff_mode is 'log_gmm'.}

\item{num_log_breaks}{Number of breaks to count frequency on the log-transformed read coverage per UMI. Only applicable when coverage_cutoff_mode is 'log_turn'.}

\item{do_plot}{Plot read coverage histogram and the threshold}
}
\description{
Filter the iTracer barcode data frame output by the extraction, based on number of supporting
reads per UMI, barcode sequences, UMI conflicts between barcodes, and UMI-barcode conflicts
between cells.
}
