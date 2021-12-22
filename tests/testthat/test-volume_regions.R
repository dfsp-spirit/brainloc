

test_that("We can compute the center of mass of coordinates.", {

    coords = matrix(c(0.0, 0.0, 0.0, 0.0, 0.0, 2.0), ncol = 3L);

    cm = brainloc:::center_of_mass(coords);

    testthat::expect_true(is.vector(cm));
    testthat::expect_equal(length(cm), 3L);

    expected = c(0.0, 0.0, 1.0);
    testthat::expect_equal(cm, expected);
})
