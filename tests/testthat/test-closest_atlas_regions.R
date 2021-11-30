
test_that("We can compute the atlas regions closest to a vertex using Euclidean distance and single linkage.", {
    fs_info = brainloc:::find.freesurferhome();
    if(! fs_info$found) {
        testthat::skip("No FreeSurfer installation found on system, but the FreeSurfer fsaverage subject is required for this test.");
    }
    sjd = file.path(fs_info$found_at, "subjects");
    sj = "fsaverage";

    options("brainloc.silent" = TRUE);

    bp = brainparc_fs(sjd, sj);
    res = vertex_closest_regions(bp, vertices=c(10, 20), hemis=c("lh", "rh"), linkage = "single", distance = "euclidean", num_regions_to_report = 3L);

    testthat::expect_true(is.data.frame(res));
    testthat::expect_equal(nrow(res), 6L); # 2 query vertices, and num_regions_to_report = 3.
    testthat::expect_equal(ncol(res), 12L); # fixed.
})


test_that("We can compute the atlas regions closest to a coordinate using Euclidean distance and single linkage.", {
    fs_info = brainloc:::find.freesurferhome();
    if(! fs_info$found) {
        testthat::skip("No FreeSurfer installation found on system, but the FreeSurfer fsaverage subject is required for this test.");
    }
    sjd = file.path(fs_info$found_at, "subjects");
    sj = "fsaverage";

    options("brainloc.silent" = TRUE);

    bp = brainparc_fs(sjd, sj);
    query_coords = get_surface(bp)$lh$vertices[1:3,];
    res = coord_closest_regions(bp, query_coords, num_regions_to_report = 5L); # We pass the coords of the first 3 vertices of the left hemisphere. Of course, you would use 'vertex_closest_regions()' if you knew the vertex index of the coordinates since it is faster.

    testthat::expect_true(is.data.frame(res));
    testthat::expect_equal(nrow(res), 15L); # 3 query vertices, and num_regions_to_report = 5.
    testthat::expect_equal(ncol(res), 12L); # fixed.
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


### Tests for all valid combinations of 'linkage' and 'distance' follow below.

test_that("We can compute the atlas regions closest to a vertex using Euclidean distance and centroid linkage.", {
    fs_info = brainloc:::find.freesurferhome();
    if(! fs_info$found) {
        testthat::skip("No FreeSurfer installation found on system, but the FreeSurfer fsaverage subject is required for this test.");
    }
    sjd = file.path(fs_info$found_at, "subjects");
    sj = "fsaverage";

    options("brainloc.silent" = TRUE);

    bp = brainparc_fs(sjd, sj);
    res = vertex_closest_regions(bp, vertices=c(10, 20), hemis=c("lh", "rh"), linkage = "centroid", distance = "euclidean", num_regions_to_report = 3L);

    testthat::expect_true(is.data.frame(res));
    testthat::expect_equal(nrow(res), 6L); # 2 query vertices, and num_regions_to_report = 3.
    testthat::expect_equal(ncol(res), 12L); # fixed.
})


test_that("We can compute the atlas regions closest to a vertex using geodesic distance and single linkage.", {
    fs_info = brainloc:::find.freesurferhome();
    if(! fs_info$found) {
        testthat::skip("No FreeSurfer installation found on system, but the FreeSurfer fsaverage subject is required for this test.");
    }
    sjd = file.path(fs_info$found_at, "subjects");
    sj = "fsaverage";

    options("brainloc.silent" = TRUE);

    bp = brainparc_fs(sjd, sj);
    res = vertex_closest_regions(bp, vertices=c(10, 20), hemis=c("lh", "rh"), linkage = "single", distance = "geodesic", num_regions_to_report = 3L);

    testthat::expect_true(is.data.frame(res));
    testthat::expect_equal(nrow(res), 6L); # 2 query vertices, and num_regions_to_report = 3.
    testthat::expect_equal(ncol(res), 12L); # fixed.
})

