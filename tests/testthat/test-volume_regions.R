

test_that("We can compute the center of mass of coordinates.", {

    coords = matrix(c(0.0, 0.0, 0.0, 0.0, 0.0, 2.0), ncol = 3L);

    cm = brainloc:::center_of_mass(coords);

    testthat::expect_true(is.vector(cm));
    testthat::expect_equal(length(cm), 3L);

    expected = c(0.0, 0.0, 1.0);
    testthat::expect_equal(cm, expected);
})


test_that("We can compute the center of mass of all atlas regions.", {

    testthat::skip_on_cran();

    if(has_fs()) {
      fsh = fs_home();
      lutfile = file.path(fsh, "FreeSurferColorLUT.txt");
      segfile = file.path(fsh, 'subjects', 'fsaverage', 'mri', 'aseg.mgz');
      named_regions = name_regions(segfile, lutfile);
      named_regions$Unknown = NULL; # Ignore an atlas region.
      sc = segmentation_centers(segfile, named_regions);
      regions_dist = dist(sc);
      dm = as.matrix(regions_dist); # convert dist object to distance matrix.

      testthat::expect_equal(ncol(dm), 43L);
      testthat::expect_equal(nrow(dm), 43L);
    }
})
