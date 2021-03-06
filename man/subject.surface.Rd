% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brainparc.R
\name{subject.surface}
\alias{subject.surface}
\title{Load a surface for a subject.}
\usage{
subject.surface(
  subjects_dir,
  subject_id,
  surface = "white",
  hemi = "both",
  force_hemilist = FALSE
)
}
\arguments{
\item{subjects_dir}{string. The FreeSurfer SUBJECTS_DIR, i.e., a directory containing the data for all your subjects, each in a subdir named after the subject identifier.}

\item{subject_id}{string. The subject identifier}

\item{surface}{string. The surface name. E.g., "white", or "pial". Used to construct the name of the surface file to be loaded.}

\item{hemi}{string, one of 'lh', 'rh', or 'both'. The hemisphere name. Used to construct the names of the surface file to be loaded. For 'both', see the information on the return value.}

\item{force_hemilist}{logical, whether to return a \code{\link{hemilist}} even if the 'hemi' parameter is not set to 'both'}
}
\value{
the `fs.surface` instance, as returned by \code{\link[freesurferformats]{read.fs.surface}}. If parameter `hemi` is set to `both`, a named list with entries `lh` and `rh` is returned, and the values of are the respective surfaces. The mesh data structure used in `fs.surface` is a *face index set*.
}
\description{
Load a brain surface mesh for a subject.
}
\concept{surface mesh functions}
\keyword{internal}
