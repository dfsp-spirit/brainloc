

#' @title Helper function for unit tests to construct a clusterinfo instance from the test data.
#'
#' @return a clusterinfo instance.
get_test_clusterinfo <- function() {
    fs_info = brainloc:::find.freesurferhome();
    if(! fs_info$found) {
        stop("No FreeSurfer installation found on system, but the FreeSurfer fsaverage subject is required for this test.");
    }
    sjd = file.path(fs_info$found_at, "subjects");
    sj = "fsaverage";


    lh_tmap_file = system.file("extdata", "lh.tmap.mgh", package = "brainloc", mustWork = TRUE);
    rh_tmap_file = system.file("extdata", "rh.tmap.mgh", package = "brainloc", mustWork = TRUE);
    lh_overlay_file = system.file("extdata", "lh.cluster.overlayID.mgh", package = "brainloc", mustWork = TRUE);
    rh_overlay_file = system.file("extdata", "rh.cluster.overlayID.mgh", package = "brainloc", mustWork = TRUE);
    clinfo = clusterinfo(lh_overlay_file, rh_overlay_file, lh_tmap_file, rh_tmap_file, template_subject = sj, subjects_dir = sjd);
    return(clinfo);
}
