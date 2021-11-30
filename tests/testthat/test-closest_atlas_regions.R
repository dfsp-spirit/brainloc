
test_that("We can compute the atlas regions closest to a vertex using Euclidean distance and single linkage.", {
    fs_info = brainloc:::find.freesurferhome();
    if(! fs_info$found) {
        testthat::skip("No FreeSurfer installation found on system, but the FreeSurfer fsaverage subject is required for this test.");
    }
    sjd = file.path(fs_info$found_at, "subjects");
    sj = "fsaverage";

    options("brainloc.silent" = TRUE);

    bp = brainparc_fs(sjd, sj);
    res = vertex_closest_regions(bp, vertices=c(10, 20), hemis=c("lh", "rh"), linkage = "single", distance = "euclidean");

    testthat::expect_equal(1L, 1L);
})



test_that("We can compute the vertex closest to the given coordinates.", {
    fs_info = brainloc:::find.freesurferhome();
    if(! fs_info$found) {
        testthat::skip("No FreeSurfer installation found on system, but the FreeSurfer fsaverage subject is required for this test.");
    }
    sjd = file.path(fs_info$found_at, "subjects");
    sj = "fsaverage";

    options("brainloc.silent" = TRUE);

    bp = brainparc_fs(sjd, sj);
    query_coords = matrix(seq(9)+0.1, ncol = 3, byrow = TRUE);
    ccv = coord_closest_vertex(query_coords, get_surface(bp));

    testthat::expect_true(is.data.frame(ccv));
    testthat::expect_equal(nrow(ccv), 3L); # for the 3 query vertices.
})
