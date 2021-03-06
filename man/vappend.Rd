% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/closest_atlas_regions.R
\name{vappend}
\alias{vappend}
\title{Append to vector or NULL.}
\usage{
vappend(vec, val)
}
\arguments{
\item{vec}{an existing vector or NULL}

\item{val}{the value to append.}
}
\value{
the value itself if vec is \code{NULL}, the vec with the appended val otherwise.
}
\description{
Append to vector or NULL.
}
\note{
There almost certainly is a better way to do this, please open an issue and let me know how if you know it.
}
\keyword{internal}
