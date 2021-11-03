

#' @title Find closest regions using Euclidean distance.
#'
#' @param brainparc a brain parcellation, see functions like \code{\link{brainparc_fs}} to get one.
#'
#' @param vertices integer vector, the query vertex indices (from the surface in the \code{brainparc}).
#'
#' @param hemis character string vector, the hemispheres for each of the \code{vertices}. Allowed entries are \code{"lh"} and \code{"rh"}, for the left and right brain hemisphere, respectively.
#'
#' @param dist_method character string, one of \code{"average"} or \code{"closest"}. Defined how the distance from a vertex to a region of vertices is computed.  \code{"average"}: Euclidean distance from query vertex to the mean of the vertex coordinates of the atlas region. \code{"closest"}: Euclidean distance from query vertex to the closest vertex of the atlas region.
#'
#' @examples
#' \dontrun{
#' bp = brainparc_fs(fsbrain::fsaverage.path(), "fsaverage", atlas="aparc");
#' vertex_closest_regions_euclid(bp, vertices=c(10, 20), hemis=c("lh", "rh"));
#' }
#' @export
vertex_closest_regions_euclid <- function(brainparc, vertices, hemis, dist_method = "average") {
    if(! (dist_method %in% c('closest', 'average'))) {
        stop("Parameter 'dist_method' must be one of c('closest', 'average').");
    }
    if(! ("brainparc" %in% class(brainparc))) {
        stop("Parameter 'brainparc' must contain a brainparc instance.");
    }
    if(length(vertices) != length(hemis)) {
        stop("Parameters 'vertices' and 'hemis' must have the same length: we need to know the hemi for each query vertex.");
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

        cat(sprintf("Handling vertex %d on hemi %s using distance method '%s'.\n", vertex_surface_idx, hemi, dist_method));

        if(dist_method == "closest") {
            vdists = freesurferformats::vertexdists.to.point(surface, vertex_coords);
            for(atlas_name in names(brainparc$annots)) {
                cat(sprintf(" -Handling atlas %s.\n", atlas_name));
                annot_min = brainparc$annots[[atlas_name]][[hemi]];
                region_names = unique(annot_min);
                vertex_region = annot_min[vertex_surface_idx];
                nr = length(region_names); # num regions

                region_idx = 0L;
                regions_closest_vertex_to_query_vertex = rep(NA, nr);
                regions_closest_distance_query_vertex = rep(NA, nr);
                for(region_name in region_names) {
                    region_idx = region_idx + 1L;
                    region_vertex_indices = which(annot_min == region_name);
                    region_vertex_dists_to_query_vertex = vdists[region_vertex_indices];
                    #sorted_region_dist_indices = sort(region_vertex_dists_to_query_vertex, index.return = TRUE)$ix;
                    closest_vertex_in_region_to_query_vertex = region_vertex_indices[which.min(region_vertex_dists_to_query_vertex)];
                    regions_closest_vertex_to_query_vertex[region_idx] = closest_vertex_in_region_to_query_vertex;
                    regions_closest_distance_query_vertex[region_idx] = euclidian.dist(vertex_coords, surface$vertices[closest_vertex_in_region_to_query_vertex, ]);
                }
                sorted_region_sort_indices = sort(regions_closest_distance_query_vertex, index.return = TRUE)$ix;
                sorted_regions = region_names[sorted_region_sort_indices];
                num_indices = min(length(sorted_regions), 5L);
                cat(sprintf("  Vertex %s on hemi %s belongs to atlas %s region '%s'. Closest region vertices are:\n", vertex_surface_idx, hemi, atlas_name, vertex_region));
                for(i in seq.int(num_indices)) {
                    cat(sprintf("  - Region #%d %s with vertex %d in distance '%f'.\n", i, region_names[sorted_region_sort_indices][i], regions_closest_vertex_to_query_vertex[sorted_region_sort_indices][i], regions_closest_distance_query_vertex[sorted_region_sort_indices][i]));
                }
            }
        } else if (dist_method == "average") {
            for(atlas_name in names(brainparc$annots)) {
                cat(sprintf(" -Handling atlas %s.\n", atlas_name));
                annot_min = brainparc$annots[[atlas_name]][[hemi]];
                region_names = unique(annot_min);
                vertex_region = annot_min[vertex_surface_idx];
                nr = length(region_names); # num regions
                region_centers_xyz = matrix(rep(NA, nr*3L), nrow = nr);
                vertex_dist_to_region_centers_xyz = rep(NA, nr);
                region_idx = 0L;
                for(region_name in region_names) {
                    region_idx = region_idx + 1L;
                    region_vertex_indices = which(annot_min == region_name);
                    region_vertex_coords = surface$vertices[region_vertex_indices, ];
                    region_center = colMeans(region_vertex_coords);
                    region_centers_xyz[region_idx, ] = region_center;
                    vertex_dist_to_region_centers_xyz[region_idx] = euclidian.dist(vertex_coords, region_center);
                    #cat(sprintf(" - Center of atlas '%s' hemi '%s' region #%d '%s' is: %f %f %f.\n", atlas_name, hemi, region_idx, region_name, region_center[1], region_center[2], region_center[3]));
                }
                sorted_region_indices = sort(vertex_dist_to_region_centers_xyz, index.return = TRUE)$ix;
                num_indices = min(length(sorted_region_indices), 5L);
                cat(sprintf("  Vertex %s on hemi %s belongs to region '%s'. Closest region centers are:\n", vertex_surface_idx, hemi, vertex_region));
                for(i in seq.int(num_indices)) {
                    cat(sprintf("  - Region #%d %s with center at (%f, %f, %f) in distance '%f'.\n", i, region_names[sorted_region_indices[i]], region_centers_xyz[sorted_region_indices[i], 1], region_centers_xyz[sorted_region_indices[i],2], region_centers_xyz[sorted_region_indices[i],3], vertex_dist_to_region_centers_xyz[sorted_region_indices[i]]));
                }
            }
        } else {
            stop("Invalid value for 'dist_method' parameter.");
        }
    }
}


#' @title Compute Euclidean distance.
#'
#' @param x1 numerical vector, coords of first point
#'
#' @param x2 numerical vector, coords of second point
#'
#' @return the Euclidean distance between x1 and x2.
#'
#' @keywords internal
euclidian.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
