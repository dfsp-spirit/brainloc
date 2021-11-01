

#' @title Find closest regions using Euclidean distance.
#'
#' @examples
#' \dontrun{
#' bp = brainparc_fs(fsbrain::fsaverage.path(), "fsaverage");
#' vertex_closest_regions_euclid(bp, c(10, 20), hemis=c("lh", "rh"));
#' }
#' @export
vertex_closest_regions_euclid <- function(brainparc, vertices, hemis, dist_method = "average") {
    if(! (dist_method %in% c('closest', 'average'))) {
        stop("Parameter 'dist_method' must be one of c('closest', 'average').");
    }

    nv = length(vertices); # nv = number of vertices.
    coords_x = rep(0.0, nv);
    coords_y = rep(0.0, nv);
    coords_z = rep(0.0, nv);
    region = rep("?", nv);


    for (vertex_local_idx in seq_along(vertices)) {
        surfaces = get_surface(brainparc);
        hemi = hemis[vertex_local_idx];
        surface = surfaces[[hemi]];
        vertex_surface_idx = vertices[vertex_local_idx];
        vertex_coords = surface$vertices[vertex_surface_idx, ]; # 1x3, xyz
        vdists = freesurferformats::vertexdists.to.point(surface, vertex_coords);


        if(dist_method == "closest") {
            stop("Sorry, dist_method 'closest' not implemented yet, try 'average' for now.");
        } else { # 'average'
            for(atlas_name in names(brainparc$annots)) {
                region_names = unique(brainparc$annots[[atlas_name]][[hemi]]);
                for(region_name in region_names) {
                    region_vertex_indices = which(brainparc$annots[[atlas_name]][[hemi]] == region_name);
                    region_vertex_coords = surface$vertices[region_vertex_indices, ];
                    region_center = colMeans(region_vertex_coords);
                    cat(sprintf("Center of atlas '%s' hemi '%s' region '%s' is: %f %f %f.\n", atlas_name, hemi, region_name, region_center[1], region_center[2], region_center[3]));
                }
            }

        }

    }

}
