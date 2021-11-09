

test_that("We can transform from MNI305 vertex to MNI152 RAS using the FreeSurfer method", {
    fs_info = brainloc:::find.freesurferhome();
    if(! fs_info$found) {
        testthat::skip("No FreeSurfer installation found on system, but tghe FreeSurfer fsaverage subject is required for this test.");
    }
    sjd = file.path(fs_info$found_at, "subjects");

    bp = brainparc_fs(sjd, "fsaverage");
    test_vertex_fsaverage_lh = 145029L;


    coord_info = coord_MNI305_info(bp$surfaces$white$lh$vertices[test_vertex_fsaverage_lh, ]);
    testthat::expect_equal(coord_info$mni152, c(-39.37023, -30.47030,  65.31317), tolerance = 0.001);
})


test_that("We can transform from MNI305 RAS to MNI152 RAS using the FreeSurfer method", {

    # This reproduces an example in the official FreeSurfer documentation, see
    # section 8 at https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems.

    fs_info = brainloc:::find.freesurferhome();
    if(! fs_info$found) {
        testthat::skip("No FreeSurfer installation found on system, but tghe FreeSurfer fsaverage subject is required for this test.");
    }
    sjd = file.path(fs_info$found_at, "subjects");

    bp = brainparc_fs(sjd, "fsaverage");
    test_RAS_MNI305 = c(10, -20, 35); # see FreeSurfer example.


    coord_info = coord_MNI305_info(test_RAS_MNI305);
    testthat::expect_equal(coord_info$mni152, c(10.695, -18.409, 36.137), tolerance = 0.001); # expected result known from FreeSurfer documentation example.
})





