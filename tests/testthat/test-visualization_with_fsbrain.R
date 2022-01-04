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

        testthat::expect_equal(1L, 1L); # Avoid skip by testthat.
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
        fsbrain::highlight.points.spheres(rbind(query_vertex_coords, region_closest_coords), color = c("red", rep("yellow", num_regions_to_report)), radius = 0.3); # draw spheres

        path_coords = matrix(rep(NA, (num_regions_to_report*2*3L)), ncol = 3L);

        local_idx = 1L;
        for(line_idx in seq.int(1L, num_regions_to_report*2, by=2L)) {
            path_coords[line_idx, ] = query_vertex_coords;
            path_coords[(line_idx+1L), ] = region_closest_coords[local_idx, ];
            local_idx = local_idx + 1L;
        }

        fsbrain::vis.paths(list(path_coords)); # draw line between spheres.

        message(sprintf("Look at brain region '%s' on the '%s' hemisphere (from top).\n", as.character(res$vertex_region_direct[1]), query_hemi));

        testthat::expect_equal(1L, 1L); # Avoid skip by testthat.
    } else {
        testthat::skip("This demo requires the optional dependency fsbrain. Please install it.");
    }
})


test_that("We can show the regions of a volume segmentation and their distances.", {
    testthat::skip_on_ci();
    testthat::skip_on_cran();

    if(requireNamespace("fsbrain", quietly = TRUE)) {
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

            reg_col = brainloc:::region_color(lutfile, unlist(unname(named_regions)));

            sjd = file.path(fsh, "subjects");
            fsbrain::vis.subject.annot(sjd, "fsaverage", atlas = "aparc", style = "semitransparent", views = "si");
            fsbrain::highlight.points.spheres(sc, color = reg_col, radius = 2.0);

            testthat::expect_equal(1L, 1L); # Avoid skip by testthat.
        }
    } else {
        testthat::skip("This demo requires the optional dependency fsbrain. Please install it.");
    }
})


test_that("We can show the connections/adjacencies between the regions of a surface parcellation.", {
    testthat::skip_on_ci();
    testthat::skip_on_cran();

    if(requireNamespace("fsbrain", quietly = TRUE)) {
        if(has_fs()) {
            fsh = fs_home();
            sj = "fsaverage";
            sjd = file.path(fsh, 'subjects');

            options("brainloc.silent" = TRUE);

            bp = brainparc_fs(sjd, sj);
            atlas = names(bp$annots); # 'aparc'

            centr = brainloc:::region_centroids(bp);
            lh_centroids = cbind(centr$lh$x, centr$lh$y, centr$lh$z);
            rh_centroids = cbind(centr$rh$x, centr$rh$y, centr$rh$z);
            region_neigh = brainparc_neighbors(bp);

            # Optional: load the full annots for the region colors.
            full_annot = list("lh"=brainloc:::subject.annot(sjd, sj, hemi="lh", atlas = atlas), "rh"=brainloc:::subject.annot(sjd, sj, hemi="rh", atlas = atlas));


            #fsbrain::vis.subject.annot(sjd, sj, atlas = atlas, style = "semitransparent", views = "si"); # draw brain
            rgl::open3d();
            rgl::bg3d("white");
            fsbrain::highlight.points.spheres(lh_centroids, color = brainloc:::region_colors(full_annot$lh, centr$lh$region), radius = 1.0); # draw spheres
            fsbrain::highlight.points.spheres(rh_centroids, color = brainloc:::region_colors(full_annot$rh, centr$rh$region), radius = 1.0); # draw spheres
            for(hemi in c("lh", "rh")) {
                neigh = region_neigh[[atlas]][[hemi]];
                for(reg1_name in colnames(neigh)) {
                    for(reg2_name in colnames(neigh)) {
                        if(nchar(reg1_name) < 1L | nchar(reg2_name) < 1L | reg1_name == "_" | reg2_name == "") {  # TODO: investigate and fix empty annot region names.
                            next;
                        }
                        if(reg1_name == reg2_name) {
                            next;
                        }
                        if(neigh[reg1_name, reg2_name] == 1L) {
                            hemi_centr = centr[[hemi]];
                            reg1_center = as.vector(hemi_centr[hemi_centr$region==reg1_name, c("x", "y", "z")]);
                            reg2_center = as.vector(hemi_centr[hemi_centr$region==reg2_name, c("x", "y", "z")]);
                            coords = list(matrix(c(reg1_center, reg2_center), ncol = 3, byrow = TRUE));
                            fsbrain::vis.paths(coords, path_color = "#777777");
                        }
                    }
                }
            }


            testthat::expect_equal(1L, 1L); # Avoid skip by testthat.
        }
    } else {
        testthat::skip("This demo requires the optional dependency fsbrain. Please install it.");
    }
})


