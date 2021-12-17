# This file does contain some examples for visualizing the clusters and their
# overlap with fsbrain. All functions require the 'fsbrain'
# package. This is more a demo than unit tests for brainloc.


test_that("We can show a coordinate and the vertex closest to it.", {
    testthat::skip_on_ci();
    testthat::skip_on_cran();
    if(requireNamespace("fsbrain", quietly = TRUE)) {
        fs_info = brainloc:::find.freesurferhome();
        if(! fs_info$found) {
            testthat::skip("No FreeSurfer installation found on system, but the FreeSurfer fsaverage subject is required for this test.");
        }
        sjd = file.path(fs_info$found_at, "subjects");
        sj = "fsaverage";

        options("brainloc.silent" = TRUE);

        bp = brainparc_fs(sjd, sj);
        query_coords = matrix(c(50,50,50), ncol = 3, byrow = TRUE);
        ccv = coord_closest_vertex(query_coords, get_surface(bp)); # find closest vertex.
        ccv_coords = get_surface_coords(bp, ccv$both_closest_vertex, ccv$both_hemi); # extract its coords.

        # Now display the source coord in red, and the closest vertex in yellow (and draw a red line between them).
        fsbrain::vis.subject.annot(sjd, sj, atlas = "aparc", views = "si"); # draw brain
        fsbrain::highlight.points.spheres(rbind(query_coords, ccv_coords), color = c("red", "yellow"), radius = 5); # draw spheres
        fsbrain::vis.paths(list(rbind(query_coords, ccv_coords))); # draw line between spheres.

        testthat::expect_equal(1L, 1L);
    } else {
        testthat::skip("This demo requires the optional dependency fsbrain. Please install it.");
    }
})


test_that("We can show a vertex and the distance to surrounding brain regions.", {
    testthat::skip_on_ci();
    testthat::skip_on_cran();
    if(requireNamespace("fsbrain", quietly = TRUE)) {
        fs_info = brainloc:::find.freesurferhome();
        if(! fs_info$found) {
            testthat::skip("No FreeSurfer installation found on system, but the FreeSurfer fsaverage subject is required for this test.");
        }
        sjd = file.path(fs_info$found_at, "subjects");
        sj = "fsaverage";

        options("brainloc.silent" = TRUE);

        bp = brainparc_fs(sjd, sj);

        query_vertex = 16564L;
        query_hemi = "lh";
        num_regions_to_report = 3L;

        res = vertex_closest_regions(bp, vertices=query_vertex, hemis=query_hemi, linkage = "single", distance = "geodesic", num_regions_to_report = num_regions_to_report);

        # Now display the query vertex and the centers of the closest regions.
        query_vertex_coords = get_surface_coords(bp, query_vertex, query_hemi);
        region_closest_coords = as.matrix(res[c("vertex_region_n_point_x","vertex_region_n_point_y","vertex_region_n_point_z")]);

        fsbrain::vis.subject.annot(sjd, sj, atlas = "aparc", views = "si"); # draw brain
        fsbrain::highlight.points.spheres(rbind(query_vertex_coords, region_closest_coords), color = c("red", rep("yellow", num_regions_to_report)), radius = 1.0); # draw spheres
        fsbrain::vis.paths(list(rbind(query_vertex_coords, region_closest_coords))); # draw line between spheres.

        message(sprintf("Look at brain region '%s' on the '%s' hemisphere (from top).\n", as.character(res$vertex_region_direct[1]), query_hemi));

        testthat::expect_equal(1L, 1L);
    } else {
        testthat::skip("This demo requires the optional dependency fsbrain. Please install it.");
    }
})

