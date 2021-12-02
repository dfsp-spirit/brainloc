

test_that("We can construct a brainparc from FreeSurfer data.", {

    sjd = get_subjects_dir(mustWork = FALSE);
    if(is.null(sjd)) {
        testthat::skip("No SUBJECTS_DIR found on system, but the FreeSurfer fsaverage subject is required for this test.");
    }

    bp = brainparc_fs(sjd, "fsaverage");
    testthat::expect_true("white" %in% names(bp$surfaces));
    testthat::expect_true("lh" %in% names(bp$surfaces$white));
    testthat::expect_true("rh" %in% names(bp$surfaces$white));

    testthat::expect_true("aparc" %in% names(bp$annots));
})
