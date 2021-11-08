

test_that("We can transform from MNI305 to Talairach using the FreeSurfer method", {
    sjd = "~/software/freesurfer/subjects";
    if(! dir.exists(sjd)) {
        testthat::skip("No test data on system.");
    }

    bp = brainparc_fs(sjd, "fsaverage");
    test_vertex_fsaverage_lh = 145029L;


    coord_info = coord_MNI305_info(bp$surfaces$white$lh$vertices[test_vertex_fsaverage_lh, ]);
    # coord_info$mni152 = -39.37023 -30.47030  65.31317
    testthat::expect_true(is.vector(coord_info)); # Just run the test to check for errors atm.
})

