% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterinfo.R
\name{clusterinfo}
\alias{clusterinfo}
\title{Create clusterinfo data structure from files or pre-loaded data.}
\usage{
clusterinfo(
  lh_overlay,
  rh_overlay,
  lh_statmap,
  rh_statmap,
  template_subject = "fsaverage",
  subjects_dir = file.path(getOption("brainloc.fs_home", default =
    Sys.getenv("FREESURFER_HOME")), "subjects")
)
}
\arguments{
\item{lh_overlay}{integer vector, the cluster overlay data: one integer per vertex of the left hemisphere. This assigns each vertex to a cluster, and all vertices of one cluster have the same number. Background is typically 1L. If a character string, the parameter will be interpreted as a file path and loaded with \code{freesurferformats::read.fs.morph}.}

\item{rh_overlay}{integer vector, just like \code{lh_overlay}, but for the right hemisphere.}

\item{lh_statmap}{double vector, the stats map. Typically a t-value map. Must have one value per vertex. If a character string, the parameter will be interpreted as file path and loaded with \code{freesurferformats::read.fs.morph}.}

\item{rh_statmap}{double vector, just like \code{lh_statmap}, but for the right hemisphere.}

\item{template_subject}{character string, the template subject name. Typically 'fsaverage' or 'fsaverage6'. Must be some MNI305 space subject.}

\item{subjects_dir}{character string, file system path to a directory containing the recon-all data for the template_subject. Used to load surfaces and annotations to identify cluster coordinates and atlas regions.}
}
\value{
named list with entries 'overlay', 'statmap', and 'metadata': a clusterinfo data structure. Each of the 'overlay' and 'statmap' keys holds a \code{\link{hemilist}} of numerical vectors.
}
\description{
Create clusterinfo data structure from files or pre-loaded data.
}
\seealso{
\code{\link{clusterinfo_from_thresholded_overlay}} If you do have a thresholded t-map instead of one t-map and one cluster overlay map.
}
