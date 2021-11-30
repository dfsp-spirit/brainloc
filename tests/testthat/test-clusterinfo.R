
test_that("We can construct a clusterinfo from separate t-maps and overlay files.", {

    fs_info = brainloc:::find.freesurferhome();
    if(! fs_info$found) {
        testthat::skip("No FreeSurfer installation found on system, but the FreeSurfer fsaverage subject is required for this test.");
    }
    sjd = file.path(fs_info$found_at, "subjects");
    sj = "fsaverage";

    options("brainloc.silent" = TRUE);
    lh_tmap_file = system.file("extdata", "lh.tmap.mgh", package = "brainloc", mustWork = TRUE);
    rh_tmap_file = system.file("extdata", "rh.tmap.mgh", package = "brainloc", mustWork = TRUE);
    lh_overlay_file = system.file("extdata", "lh.cluster.overlayID.mgh", package = "brainloc", mustWork = TRUE);
    rh_overlay_file = system.file("extdata", "rh.cluster.overlayID.mgh", package = "brainloc", mustWork = TRUE);
    clinfo = clusterinfo(lh_overlay_file, rh_overlay_file, lh_tmap_file, rh_tmap_file, template_subject = sj, subjects_dir = sjd);
    testthat::expect_true(is.clusterinfo(clinfo));
    testthat::expect_equal(length(get_clusters(clinfo)), 2L); # We know there are 2 clusters in the data.
})


test_that("We can construct a clusterinfo from thresholded t-maps.", {

    fs_info = brainloc:::find.freesurferhome();
    if(! fs_info$found) {
        testthat::skip("No FreeSurfer installation found on system, but the FreeSurfer fsaverage subject is required for this test.");
    }
    sjd = file.path(fs_info$found_at, "subjects");
    sj = "fsaverage";

    options("brainloc.silent" = TRUE);
    lh_tmap_file = system.file("extdata", "lh.tmap.mgh", package = "brainloc", mustWork = TRUE);
    rh_tmap_file = system.file("extdata", "rh.tmap.mgh", package = "brainloc", mustWork = TRUE);
    lh_overlay_file = system.file("extdata", "lh.cluster.overlayID.mgh", package = "brainloc", mustWork = TRUE);
    rh_overlay_file = system.file("extdata", "rh.cluster.overlayID.mgh", package = "brainloc", mustWork = TRUE);

    # We construct the tresholded t-map from full data for this example, which one would not do for real data of course.
    # It is way simpler and faster to call 'clusterinfo()' directly if you have both the maps and an overlay.
    lh_tmap = freesurferformats::read.fs.morph(lh_tmap_file);
    rh_tmap = freesurferformats::read.fs.morph(rh_tmap_file);
    lh_overlay = as.integer(freesurferformats::read.fs.morph(lh_overlay_file));
    rh_overlay = as.integer(freesurferformats::read.fs.morph(rh_overlay_file));

    # Construct thresholded map: set the t-map values of all vertices which are not in any cluster to 0.
    lh_threshmap = lh_tmap;
    rh_threshmap = rh_tmap;
    lh_threshmap[lh_overlay == 0] = 0.0;
    rh_threshmap[rh_overlay == 0] = 0.0;

    clinfo = clusterinfo_from_thresholded_overlay(lh_threshmap, rh_threshmap, template_subject = sj, subjects_dir = sjd);
    testthat::expect_true(is.clusterinfo(clinfo));
    testthat::expect_equal(length(get_clusters(clinfo)), 2L); # We know there are 2 clusters in the data.
})
