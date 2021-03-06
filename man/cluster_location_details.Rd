% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster_loc.R
\name{cluster_location_details}
\alias{cluster_location_details}
\title{Get details on cluster locations.}
\usage{
cluster_location_details(
  clusterinfo,
  silent = getOption("brainloc.silent", default = FALSE),
  ...
)
}
\arguments{
\item{clusterinfo}{a clusterinfo instance, see the \code{clusterinfo} function on how to get one. Must contain a valid \code{brainparc} in field 'brainparc'.}

\item{silent}{logical, whether to suppress console messages.}

\item{...}{passed on to \code{\link{cluster_extrema}}.}
}
\value{
a new version of the input data.frame, with additional columns appended.
}
\description{
Get details on cluster locations.
}
