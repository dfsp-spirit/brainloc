

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
#'
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

    num_regions_to_report = 3L;

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
                num_indices = min(length(sorted_regions), num_regions_to_report);
                cat(sprintf("  Vertex %s on hemi %s at (%f %f %f) belongs to atlas %s region '%s'. Closest region vertices are:\n", vertex_surface_idx, hemi, vertex_coords[1], vertex_coords[2], vertex_coords[3], atlas_name, vertex_region));
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
                num_indices = min(length(sorted_region_indices), num_regions_to_report);
                cat(sprintf("  Vertex %s on hemi %s at (%f %f %f) belongs to region '%s'. Closest region centers are:\n", vertex_surface_idx, hemi, vertex_coords[1], vertex_coords[2], vertex_coords[3], vertex_region));
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



#' @title Find closest surface vertex to a point, using Euclidean distance.
#'
#' @param coordinate \code{nx3} numerical matrix or vector of length 3, the query point coordinates.
#'
#' @param surfaces hemilist of \code{fs.surface} instances
#'
#' @return a data.frame with columns named 'query_x', 'query_y', 'query_z', 'lh_closest_vertex', 'lh_distance', 'rh_closest_vertex', 'rh_distance', 'both_closest_vertex', 'both_distance', 'both_hemi'.
#'
#' #' @examples
#' \dontrun{
#' bp = brainparc_fs(fsbrain::fsaverage.path(), "fsaverage", atlas="aparc");
#' query_coords = matrix(seq.int(9), ncol = 3, byrow = TRUE);
#' coord_closest_vertex(query_coords, bp$surfaces$white);
#' }
#'
#' @export
coord_closest_vertex <- function(coordinate, surfaces) {
    if(is.vector(coordinate) & length(coordinate) == 3L) {
        coordinate = matrix(coordinate, ncol = 3, nrow = 1, byrow = TRUE);
    }
    if(! is.matrix(coordinate)) {
        stop("Parameter 'coordinate' must be a vector of length 3 or an nx3 matrix.");
    }

    num_coords = nrow(coordinate);
    lh_closest_vertex = rep(NA, num_coords);
    rh_closest_vertex = rep(NA, num_coords);
    both_closest_vertex = rep(NA, num_coords);
    lh_distance = rep(NA, num_coords);
    rh_distance = rep(NA, num_coords);
    both_distance = rep(NA, num_coords);
    both_hemi = rep(NA, num_coords);

    for (row_idx in seq.int(num_coords)) {
        has_surf = FALSE;
        if(freesurferformats::is.fs.surface(surfaces$lh)) {
            lh_vd = freesurferformats::vertexdists.to.point(surfaces$lh, coordinate[row_idx, ]);
            lh_closest_vertex[row_idx] = which.min(lh_vd);
            lh_distance[row_idx] = lh_vd[lh_closest_vertex[row_idx]];
            both_closest_vertex[row_idx] = lh_closest_vertex[row_idx]; # for now, may change below.
            both_distance[row_idx] = lh_distance[row_idx];             # for now, may change below.
            both_hemi[row_idx] = "lh";
            has_surf = TRUE;
        }
        if(freesurferformats::is.fs.surface(surfaces$rh)) {
            rh_vd = freesurferformats::vertexdists.to.point(surfaces$rh, coordinate[row_idx, ]);
            rh_closest_vertex[row_idx] = which.min(rh_vd);
            rh_distance[row_idx] = rh_vd[rh_closest_vertex[row_idx]];
            if(has_surf) { # whether there was a left surface. in this case we need to compare the values and pick the smaller dist.
                if(lh_distance[row_idx] < rh_distance[row_idx]) {
                    both_closest_vertex[row_idx] = lh_closest_vertex[row_idx];
                    both_distance[row_idx] = lh_distance[row_idx];
                    both_hemi[row_idx] = "lh";
                } else {
                    both_closest_vertex[row_idx] = rh_closest_vertex[row_idx];
                    both_distance[row_idx] = rh_distance[row_idx];
                    both_hemi[row_idx] = "rh";
                }
            } else { # Only the rh surface was given.
                both_closest_vertex[row_idx] = rh_closest_vertex[row_idx];
                both_distance[row_idx] = rh_distance[row_idx];
                both_hemi[row_idx] = "rh";
            }
            has_surf = TRUE;
        }
        if(! has_surf) {
            stop("The hemilist in parameter 'surfaces' must contain at least one fs.surface instance in keys 'lh' or 'rh'.");
        }
    }
    return(data.frame("query_x"=coordinate[,1], "query_y"=coordinate[,2], "query_z"=coordinate[,3], "lh_closest_vertex"=lh_closest_vertex, "lh_distance"=lh_distance, "rh_closest_vertex"=rh_closest_vertex, "rh_distance"=rh_distance, "both_closest_vertex"=both_closest_vertex, "both_distance"=both_distance, "both_hemi"=both_hemi));
}

