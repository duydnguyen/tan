% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{plotCoverage}
\alias{plotCoverage}
\title{create an enrichment profile for a set of samples at a given peak}
\usage{
plotCoverage(coverage, peak_id, sample_ids = NA, title = "", size = 0.5,
  showLegend = TRUE)
}
\arguments{
\item{coverage}{: A list of coverage whose elements contain reads from given
regions.}

\item{peak_id}{: integer indicating the index of the peak to be plotted}

\item{sample_ids}{: labels of the considered samples for enrichment plot. If
NA, all samples will be used, and the mean profiles for each conditions
are created.}

\item{title}{: title of the plot}

\item{size}{: thickness of lines}

\item{showLegend}{: if TRUE legend is shown}
}
\value{
A ggplot2 plot
}
\description{
create an enrichment profile for a set of samples at a given peak
}
