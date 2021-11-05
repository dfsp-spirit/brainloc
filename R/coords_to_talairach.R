

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
