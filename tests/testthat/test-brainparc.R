

test_that("We can construct a brainparc from FreeSurfer data.", {

    fs_info = brainloc:::find.freesurferhome();
    if(! fs_info$found) {
        testthat::skip("No FreeSurfer installation found on system, but tghe FreeSurfer fsaverage subject is required for this test.");
    }
    sjd = file.path(fs_info$found_at, "subjects");

    bp = brainparc_fs(sjd, "fsaverage");
    testthat::expect_true("white" %in% names(bp$surfaces));
    testthat::expect_true("lh" %in% names(bp$surfaces$white));
    testthat::expect_true("rh" %in% names(bp$surfaces$white));

    testthat::expect_true("aparc" %in% names(bp$annots));
})
