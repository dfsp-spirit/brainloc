% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_loc.R
\name{cluster_peaks}
\alias{cluster_peaks}
\title{Compute cluster peaks.}
\usage{
cluster_peaks(
  clusterinfo,
  type = "extreme",
  silent = getOption("brainloc.silent", default = FALSE),
  ...
)
}
\arguments{
\item{clusterinfo}{a clusterinfo instance, see the \code{clusterinfo} function to get one.}

\item{type}{character string, one of 'extreme', 'min' or 'max'. Which cluster value to report. The default of 'extreme' reports the absolutely larger one of the min and the max value and is for the typical use case of \code{t}-value maps.}

\item{silent}{whether to suppress printing messages to stdout.}

\item{...}{passed on to \code{\link{clusteroverlay_to_annot}}.}
}
\value{
a \code{data.frame} with cluster peak information. The column names should be self-explanatory.
}
\description{
A peak is a position in a cluster (a vertex) which is assigned the maximum statsmap value over its neighborhood. The neighborhood definition used by this function is the 1-ring neighborhood on the mesh.
}
