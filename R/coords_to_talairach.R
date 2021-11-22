

#' @title Get info on MNI305 surface coordinate in various coordinate systems.
#'
#' @param coords_mni305 nx3 matrix of coordinates, typically from fsaverage surface vertices.
#'
#' @param ... Passed on to \code{\link{coord_MNI305_to_MNI152}}.
#'
#' @return named list with coordinate information (matrices)
#'
#' @export
coord_MNI305_info <- function(coords_mni305, ...) {
  coords_mni152 = coord_MNI305_to_MNI152(coords_mni305, ...);
  return(list("mni305" = coords_mni305, "mni152" = coords_mni152, "talairach" = coord_MNI152_to_talairach(coords_mni152)));
}


#' @title Transform FreeSurfer surface space coordinates into MNI Talairach space using talairach.xfm file.
#'
#' @description This uses the affine transform in \code{<subject>/mri/transforms/talairach.xfm}, be sure to know about the limitations of this approach. See the FreeSurfer documentation on talairach for details. The coordinates produced by this approach are referred to as 'MNI Talairach' coordinates in the FreeSurfer documentation, and are not identical to Talairach space.
#'
#' @param subjects_dir character string, the path to the recon-all SUBJECTS_DIR containing the data of all subjects of your study.
#'
#' @param subject_id character string, the subject identifier (i.e., the sub directory name under the \code{subjects_dir}). This can be a template like fsaverage, but it also works for other subjects.
#'
#' @param surface_coords nx3 numerical matrix of surface coordinates (vertex positions in FreeSurfer surface space).
#'
#' @return nx3 numerical matrix of MNI Talairach space coordinates.
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
#' @keywords internal
coord_fssurface_to_fstalairach <- function(subjects_dir, subject_id, surface_coords) {
    talairach_file = file.path(subjects_dir, subject_id, "mri", "transforms", "talairach.xfm");
    if(! file.exists(talairach_file)) {
      stop(sprintf("Could not read talairach file for subject '%s' at '%s'.\n", subject_id, talairach_file));
    }
    tal_coords = freesurferformats::ras.to.talairachras(surface_coords, talairach_file);
    return(tal_coords);
}


#' @title Transform MNI305 coords (FreeSurfer fsaverage surface) to MNI152 coordinates.
#'
#' @param vertex_coords \code{nx3} numerical matrix of MNI305 surface RAS coordinates, typically from fsaverage surface vertices.
#'
#' @param method character string, the method to use to map from MNI305 to MNI152 along the way. One of "best_available", "regfusionr", and "linear".
#'
#' @param surface optional character string or hemilist of surfaces, the surface to use to find a surface vertex close to the given query coordinates. Only used if 'method' results in \code{regfusionr} being used. Passed on to \code{regfusionrregfusionr::mni305_coords_to_mni152_coords}.
#'
#' @param fs_home optional character string, the path to the FREESURFER_HOME directory from which to load the surfaces from the 'surface' parameter. Only used if 'method' results in \code{regfusionr} being used. Passed on to \code{regfusionrregfusionr::mni305_coords_to_mni152_coords}.
#'
#' @return \code{nx3} numerical matrix of MNI152 coords.
#'
#' @note One can check that the results are okay by clicking a vertex in FreeView, using the displayed MNI305 coordinates as input to this function, and looking up the reported MNI152 coordinates at \code{https://bioimagesuiteweb.github.io/webapp/mni2tal.html}.
#'
#' @examples
#' \dontrun{
#' # Get MNI152 coordinates for first 3 fsaverage lh vertices:
#' surf = freesurferformats::read.fs.surface("/opt/freesurfer/subjects/fsaverage/surf/lh.white");
#' coord_MNI305_to_MNI152(surf$vertices[1:3, ]);
#' }
#'
#' @export
coord_MNI305_to_MNI152 <- function(vertex_coords, method = getOption("brainloc.method_MNI305_to_from_MNI152", default="linear"), surface = "orig", fs_home = getOption("brainloc.fs_home", default = Sys.getenv("FREESURFER_HOME"))) {
  if(! (method %in% c("best_available", "regfusionr", "linear"))) {
    stop("Parameter 'method' must be one of c('best_available', 'regfusionr', 'linear').");
  }
  if(method %in% c("best_available", "regfusionr")) {
    if(requireNamespace("regfusionr", quietly = TRUE)) {
      if(! is.list(surface)) {
        surface_name = surface;
        fsaverage_path = file.path(fs_home, 'subjects', 'fsaverage');
        if(! dir.exists(fsaverage_path)) {
          stop(sprintf("The fsaverage data for regfusionr cannot be loaded: directory '%s' does not exist or cannot be read.\n", fsaverage_path));
        } else {
          lh_surf = freesurferformats::read.fs.surface(file.path(fsaverage_path, "surf", sprintf("lh.%s", surface_name)));
          rh_surf = freesurferformats::read.fs.surface(file.path(fsaverage_path, "surf", sprintf("rh.%s", surface_name)));
          surface = list("lh" = lh_surf, "rh" = rh_surf);
        }
      }
      return(regfusionr::mni305_coords_to_mni152_coords(vertex_coords, surface = surface, fs_home = fs_home, simplify = TRUE));
    } else {
      if(method == "regfusionr") {
        stop("Parameter 'method' forces regfusionr but package not available. Please install the regfusionr package from https://github.com/dfsp-spirit/regfusionr or change the 'method' parameter.");
      }
    }
  }
  return(freesurferformats::doapply.transform.mtx(vertex_coords, mni152reg_mtx()));
}


#' @title Get fsaverage (MNI305) to MNI152 transformation matrix.
#'
#' @description This returns the 4x4 matrix from the FreeSurfer Coordinate Systems documentation to transform between MNI305 and MNI152.
#'
#' @keywords internal
mni152reg_mtx <- function() {
  return(matrix(c(0.9975, -0.0073, 0.0176, -0.0429, 0.0146, 1.0009, -0.0024, 1.5496, -0.0130, -0.0093, 0.9971, 1.1840, 0, 0, 0, 1), ncol = 4, byrow = TRUE));
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
#' @examples
#'     mni_coords = matrix(c(10, 12, 14), nrow = 1, ncol = 3, byrow = TRUE);
#'     coord_MNI152_to_talairach(mni_coords);
#'
#' @export
coord_MNI152_to_talairach <- function(mni152_coords) {

  was_vector = FALSE;
  if(! is.matrix(mni152_coords)) {
    was_vector = TRUE;
    if(is.vector(mni152_coords) & length(mni152_coords) == 3L) {
      mni152_coords = matrix(mni152_coords, ncol = 3, byrow = TRUE);
    } else {
      stop("Parameter 'mni152_coords' must be a matrix.");
    }
  }

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

  tal_coords = matrix(rep(NA, length(mni152_coords)), ncol = 3);
  for(row_idx in seq.int(nrow(mni152_coords))) {
    if(mni152_coords[row_idx, 3] < 0) { # whether below AZ
      tal_coords[row_idx, ] = freesurferformats::doapply.transform.mtx(mni152_coords[row_idx, ], (mtx_MNI152toTal_rotn %*% mtx_MNI152toTal_downZ));
    } else {
      tal_coords[row_idx, ] = freesurferformats::doapply.transform.mtx(mni152_coords[row_idx, ], (mtx_MNI152toTal_rotn %*% mtx_MNI152toTal_upZ));
    }
  }
  if(was_vector) {
    tal_coords = as.double(tal_coords);
  }
  return(tal_coords);
}

