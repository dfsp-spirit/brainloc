

#' @title Transform FreeSurfer surface space coordinates into Talairach space.
#'
#' @description This uses the affine transform in \code{<subject>/mri/transforms/talairach.xfm}, be sure to know about the limitations of this approach. See the FreeSurfer documentation on talairach for details.
#'
#' @param subjects_dir character string, the path to the recon-all SUBJECTS_DIR containing the data of all subjects of your study.
#'
#' @param subject_id character string, the subject identifier (i.e., the sub directory name under the \code{subjects_dir}).
#'
#' @param surface_coords nx3 numerical matrix of surface coordinates (vertex positions in FreeSurfer surface space).
#'
#' @return nx3 numerical matrix of Talairach space coordinates.
#'
#' @examples
#' \dontrun{
#' # Will lead to identity transform for fsaverage, so no change in coords.
#' sjd = fsbrain::fsaverage.path();
#' sj = "fsaverage";
#' surf = fsbrain::subject.surface(sjd, sj, hemi="lh", surface="orig");
#' tal = coord_fssurface_to_talairach(sjd, sj, surf$vertices);
#' }
#'
#' @export
coord_fssurface_to_talairach <- function(subjects_dir, subject_id, surface_coords) {
    talairach_file = file.path(subjects_dir, subject_id, "mri", "transforms", "talairach.xfm");
    if(! file.exists(talairach_file)) {
      stop(sprintf("Could not read talairach file for subject '%s' at '%s'.\n", subject_id, talairach_file));
    }
    tal_coords = freesurferformats::ras.to.talairachras(surface_coords, talairach_file);
    return(tal_coords);
}


#' @title Transform MNI305 coords (FreeSurfer fsaverage surface) to MNI152 coordinates.
#'
#' @param vertex_coords nx3 matrix of coordinates, e.g., typically from fsaverage surface vertices.
#'
#' @return nx3 numerical matrix if MNI152 coords.
#'
#' @keywords internal
coord_fsaverage_to_MNI152 <- function(vertex_coords) {
  return(freesurferformats::doapply.transform.mtx(vertex_coords, freesurferformats::mni152reg()));
}


#' @title Transform MNI152 coords to Talairach space using Matthew Brett's approach.
#'
#' @description See \code{http://brainmap.org/training/BrettTransform.html}.
#'
#' @param mni152_coords nx3 numerical matrix of RAS coordinates in MNI152 space.
#'
#' @return nx3 numerical matrix of Talairach coordinates.
#'
#' @note This is published under the GPL license. See \code{https://github.com/sccn/dipfit/blob/master/mni2tal_matrix.m} and \code{https://github.com/sccn/dipfit/blob/master/mni2tal.m} for a Matlab implementation of the method. All credits go to Matthew Brett.
#'
#' @keywords internal
coord_MNI152_to_talairach <- function(mni152_coords) {
  mtx_MNI152toTal_rotn  = matrix(c(1, 0, 0, 0,
                                 0, 0.9988, 0.0500, 0,
                                 0, -0.0500, 0.9988, 0,
                                 0, 0, 0, 1.0000), nrow = 4, ncol = 4, byrow = TRUE);

  mtx_MNI152toTal_upZ = matrix(c(0.9900, 0, 0, 0,
                0, 0.9700, 0, 0,
                0, 0, 0.9200, 0,
                0, 0, 0, 1.0000), nrow = 4, ncol = 4, byrow = TRUE);

  mtx_MNI152toTal_downZ = matrix(c(0.9900, 0, 0, 0,
                0,    0.9700,         0,         0,
                0,         0,    0.8400,         0,
                0,         0,         0,    1.0000), nrow = 4, ncol = 4, byrow = TRUE);
}

