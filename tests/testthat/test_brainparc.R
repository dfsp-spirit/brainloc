

test_that("We can construct a brainparc from FreeSurfer data.", {

    sjd = "~/software/freesurfer/subjects";
    if(! dir.exists(sjd)) {
        testthat::skip("No test data on system.");
    }

    bp = brainparc_fs(sjd, "fsaverage");
    testthat::expect_true("white" %in% names(bp$surfaces));
    testthat::expect_true("lh" %in% names(bp$surfaces$white));
    testthat::expect_true("rh" %in% names(bp$surfaces$white));

    testthat::expect_true("aparc" %in% names(bp$annots$white));
})
