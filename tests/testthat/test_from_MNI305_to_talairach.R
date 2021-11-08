

test_that("We can transform from MNI305 to Talairach using the FreeSurfer method", {
    sjd = "~/software/freesurfer/subjects";
    if(! dir.exists(sjd)) {
        testthat::skip("No test data on system.");
    }

    bp = brainparc_fs(sjd, "fsaverage");
    coord_MNI305_info()
})
