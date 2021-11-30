
test_that("We can compute cluster location details from a clusterinfo instance.", {

    options("brainloc.silent" = TRUE);
    clinfo = get_test_clusterinfo();
    testthat::expect_true(is.clusterinfo(clinfo));

    extrema_details = cluster_location_details(clinfo);

    testthat::expect_true(is.data.frame(extrema_details));
    testthat::expect_equal(nrow(extrema_details), 2L); # 1 row per cluster.
})


test_that("We can compute cluster peaks from a clusterinfo instance.", {

    options("brainloc.silent" = TRUE);
    clinfo = get_test_clusterinfo();
    testthat::expect_true(is.clusterinfo(clinfo));

    peaks = cluster_peaks(clinfo);

    testthat::expect_true(is.data.frame(peaks));
    testthat::expect_equal(nrow(peaks), 769L);
})

