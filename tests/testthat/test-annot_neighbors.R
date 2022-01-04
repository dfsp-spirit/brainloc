

test_that("We can compute the neighborhood of brain surface regions.", {

    sjd = get_subjects_dir(mustWork = FALSE);
    if(is.null(sjd)) {
        testthat::skip("No SUBJECTS_DIR found on system, but the FreeSurfer fsaverage subject is required for this test.");
    }

    bp = brainparc_fs(sjd, "fsaverage");

    annot_min = bp$annots$aparc$lh;
    surface = bp$surfaces$white$lh

    region_neigh = brainloc:::annot_neighbors(annot_min, surface);

    testthat::expect_true(is.matrix(region_neigh));

    frontal_pole_neighbors = colnames(region_neigh)[region_neigh["frontalpole",]==1L];
    testthat::expect_true("frontalpole" %in% frontal_pole_neighbors);
    testthat::expect_equal(length(frontal_pole_neighbors), 5L);
})
