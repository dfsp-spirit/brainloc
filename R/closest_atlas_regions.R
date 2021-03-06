# Functions to compute the atlas regions closest to a point/vertex coordinate on the surface.

#' @title Find regions for vertices in all brainparc atlases.
#'
#' @inheritParams vertex_closest_regions
#'
#' @return a data.frame, the column names should be obvious.
#'
#' @examples
#' \dontrun{
#' bp = brainparc_fs(get_subjects_dir(), "fsaverage", atlas="aparc");
#' vertex_region(bp, vertices=c(10, 20), hemis=c("lh", "rh"));
#' }
#'
#' @seealso \code{\link{vertex_closest_regions}} to compute distance to all regions.
#'
#' @export
vertex_regions <- function(brainparc, vertices, hemis) {
    if(! is.brainparc(brainparc)) {
        stop("Parameter 'brainparc' must contain a brainparc instance.");
    }
    nv = length(vertices);
    if(nv < 1L) {
        stop("Parameter 'vertices' must contain at least one vertex.");
    }
    if(nv != length(hemis)) {
        if(length(hemis) == 1L) {
            hemis = rep(hemis, nv); # Assume all the same hemi.
        } else {
            stop("Parameters 'vertices' and 'hemis' must have the same length: we need to know the hemi for each query vertex.");
        }
    }

    atlas = rep("?", nv);
    vertex_region = rep("?", nv);

    for(vidx in seq_along(vertices)) {
        vertex_surface_idx = vertices[vidx];
        hemi = hemis[vidx];
        for(atlas_name in names(brainparc$annots)) {
            annot_min = brainparc$annots[[atlas_name]][[hemi]];
            atlas[vidx] = atlas_name;
            vertex_region[vidx] = annot_min[vertex_surface_idx];
        }
    }
    return(data.frame("vertex"=vertices, "hemi"=hemis, "atlas"=atlas, "vertex_region"=vertex_region, stringsAsFactors = FALSE));
}


#' @title Find closest regions to vertex using Euclidean or geodesic distance.
#'
#' @description Finds the closest atlas regions according to the brain surface parcellation 'brainparc' for the given query vertices.
#'
#' @param brainparc a brain parcellation, see functions like \code{\link{brainparc_fs}} to get one.
#'
#' @param vertices integer vector, the query vertex indices (from the surface in the \code{brainparc}).
#'
#' @param hemis character string vector, the hemispheres for each of the \code{vertices}. Allowed entries are \code{"lh"} and \code{"rh"}, for the left and right brain hemisphere, respectively. Must have same length as the 'vertices' vector (or exactly length 1, in which case we assume that this is the hemi for ALL query vertices).
#'
#' @param linkage character string, one of \code{"single"} or \code{"centroid"}. Defines how the distance from a vertex to a region of vertices is computed. \code{"single"}: Euclidean distance from query vertex to the closest vertex of the atlas region. \code{"centroid"}: Euclidean distance from query vertex to the mean of the vertex coordinates of the atlas region.
#'
#' @param distance character string, one of \code{"euclidean"} or \code{"geodesic"}. The latter is only supported with \code{linkage = 'single'}.
#'
#' @param silent logical, whether to suppress console messages.
#'
#' @param num_regions_to_report scalar integer, the number of closest regions to report (will be reported in increasing order, by distance). Pass \code{NULL} to get all regions.
#'
#' @return a data.frame, the column names should be obvious.
#'
#' @examples
#' \dontrun{
#' bp = brainparc_fs(get_subjects_dir(), "fsaverage", atlas="aparc");
#' vertex_closest_regions(bp, vertices=c(10, 20), hemis=c("lh", "rh"));
#' }
#'
#' @seealso \code{\link{coord_closest_regions}} if you have a coordinate (on or near the surface) instead of a vertex.
#'
#' @export
vertex_closest_regions <- function(brainparc, vertices, hemis, linkage = "single", distance = "euclidean", silent = getOption("brainloc.silent", default = FALSE), num_regions_to_report = 5L) {
    if(! (linkage %in% c('single', 'centroid'))) {
        stop("Parameter 'linkage' must be one of c('single', 'centroid').");
    }
    if(! (distance %in% c('euclidean', 'geodesic'))) {
        stop("Parameter 'distance' must be one of c('euclidean', 'geodesic').");
    }
    if(distance == "geodesic" & linkage == "centroid") {
        stop("The distance type 'geodesic' is only supported with linkage = 'single'.");
    }
    if(! is.brainparc(brainparc)) {
        stop("Parameter 'brainparc' must contain a brainparc instance.");
    }
    if(length(vertices) != length(hemis)) {
        if(length(hemis) == 1L) {
            hemis = rep(hemis, length(vertices)); # Assume all the same hemi.
        } else {
            stop("Parameters 'vertices' and 'hemis' must have the same length: we need to know the hemi for each query vertex.");
        }
    }

    # Will be turned into data.frame later.
    res = list("vertex"=NULL, "hemi"=NULL, "linkage"=NULL, "distance"=NULL, "atlas"=NULL, "vertex_region_direct"=NULL, "vertex_region_nth"=NULL, "vertex_region_n_name"=NULL, "vertex_region_n_distance"=NULL, "vertex_region_n_point_x"=NULL, "vertex_region_n_point_y"=NULL, "vertex_region_n_point_z"=NULL);


    for (vertex_local_idx in seq_along(vertices)) {
        surfaces = get_surface(brainparc);
        hemi = hemis[vertex_local_idx];
        surface = surfaces[[hemi]];
        vertex_surface_idx = vertices[vertex_local_idx];
        vertex_coords = surface$vertices[vertex_surface_idx, ]; # 1x3, xyz

        if(! silent) {
            cat(sprintf("Handling vertex %d on hemi %s using %s distance and %s linkage.\n", vertex_surface_idx, hemi, distance, linkage));
        }

        if(linkage == "single") {
            vdists = NULL;
            if(distance == "geodesic") {
                tmesh = fs.surface.to.tmesh3d(surface);
                vdists = Rvcg::vcgDijkstra(tmesh, vertex_surface_idx);
            } else if (distance == "euclidean") {
                vdists = freesurferformats::vertexdists.to.point(surface, vertex_coords);
            } else {
                stop("Invalid 'distance' parameter.")
            }
            for(atlas_name in names(brainparc$annots)) {
                if(! silent) {
                    cat(sprintf(" -Handling atlas %s.\n", atlas_name));
                }
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
                    local_regions_closest_vertex_to_query_vertex = which.min(region_vertex_dists_to_query_vertex);
                    closest_vertex_in_region_to_query_vertex = region_vertex_indices[local_regions_closest_vertex_to_query_vertex];
                    regions_closest_vertex_to_query_vertex[region_idx] = closest_vertex_in_region_to_query_vertex;
                    regions_closest_distance_query_vertex[region_idx] = min(region_vertex_dists_to_query_vertex);
                }
                sorted_region_sort_indices = sort(regions_closest_distance_query_vertex, index.return = TRUE)$ix;
                sorted_regions = region_names[sorted_region_sort_indices];
                num_reg = num_regions_to_report;
                if(is.null(num_reg)) {
                    num_reg = length(sorted_regions);
                }
                num_indices = min(length(sorted_regions), num_reg);
                if(! silent) {
                    cat(sprintf("  Vertex %s on hemi %s at (%f %f %f) belongs to atlas %s region '%s'. Closest region vertices with %s distance are:\n", vertex_surface_idx, hemi, vertex_coords[1], vertex_coords[2], vertex_coords[3], atlas_name, vertex_region, distance));
                }
                for(i in seq.int(num_indices)) {
                    if(! silent) {
                        cat(sprintf("  - Region #%d %s with vertex %d in %s distance '%f'.\n", i, region_names[sorted_region_sort_indices][i], regions_closest_vertex_to_query_vertex[sorted_region_sort_indices][i], distance, regions_closest_distance_query_vertex[sorted_region_sort_indices][i]));
                    }

                    # Collect results.
                    res$vertex = vappend(res$vertex, vertex_surface_idx);
                    res$hemi = vappend(res$hemi, hemi);
                    res$linkage = vappend(res$linkage, linkage);
                    res$distance = vappend(res$distance, distance);
                    res$atlas = vappend(res$atlas, atlas_name);
                    res$vertex_region_direct = vappend(res$vertex_region_direct, vertex_region);
                    res$vertex_region_nth = vappend(res$vertex_region_nth, i);
                    res$vertex_region_n_name = vappend(res$vertex_region_n_name, region_names[sorted_region_sort_indices][i]);
                    res$vertex_region_n_distance = vappend(res$vertex_region_n_distance, regions_closest_distance_query_vertex[sorted_region_sort_indices][i]);

                    region_measure_vertex = regions_closest_vertex_to_query_vertex[sorted_region_sort_indices][i];

                    res$vertex_region_n_point_x = vappend(res$vertex_region_n_point_x, surface$vertices[region_measure_vertex, 1]);
                    res$vertex_region_n_point_y = vappend(res$vertex_region_n_point_y, surface$vertices[region_measure_vertex, 2]);
                    res$vertex_region_n_point_z = vappend(res$vertex_region_n_point_z, surface$vertices[region_measure_vertex, 3]);
                }
            }
        } else if (linkage == "centroid") {
            for(atlas_name in names(brainparc$annots)) {
                if(! silent) {
                    cat(sprintf(" -Handling atlas %s.\n", atlas_name));
                }
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
                }
                sorted_region_indices = sort(vertex_dist_to_region_centers_xyz, index.return = TRUE)$ix;
                num_reg = num_regions_to_report;
                if(is.null(num_reg)) {
                    num_reg = length(sorted_region_indices);
                }
                num_indices = min(length(sorted_region_indices), num_reg);
                if(! silent) {
                    cat(sprintf("  Vertex %s on hemi %s at (%f %f %f) belongs to region '%s'. Closest region centers are:\n", vertex_surface_idx, hemi, vertex_coords[1], vertex_coords[2], vertex_coords[3], vertex_region));
                }
                for(i in seq.int(num_indices)) {
                    if(! silent) {
                        cat(sprintf("  - Region #%d %s with center at (%f, %f, %f) in distance '%f'.\n", i, region_names[sorted_region_indices[i]], region_centers_xyz[sorted_region_indices[i], 1], region_centers_xyz[sorted_region_indices[i],2], region_centers_xyz[sorted_region_indices[i],3], vertex_dist_to_region_centers_xyz[sorted_region_indices[i]]));
                    }

                    # Collect results.
                    res$vertex = vappend(res$vertex, vertex_surface_idx);
                    res$hemi = vappend(res$hemi, hemi);
                    res$linkage = vappend(res$linkage, linkage);
                    res$distance = vappend(res$distance, distance);
                    res$atlas = vappend(res$atlas, atlas_name);
                    res$vertex_region_direct = vappend(res$vertex_region_direct, vertex_region);
                    res$vertex_region_nth = vappend(res$vertex_region_nth, i);
                    res$vertex_region_n_name = vappend(res$vertex_region_n_name, region_names[sorted_region_indices[i]]);
                    res$vertex_region_n_distance = vappend(res$vertex_region_n_distance, vertex_dist_to_region_centers_xyz[sorted_region_indices[i]]);
                    res$vertex_region_n_point_x = vappend(res$vertex_region_n_point_x, region_centers_xyz[sorted_region_indices[i], 1]);
                    res$vertex_region_n_point_y = vappend(res$vertex_region_n_point_y, region_centers_xyz[sorted_region_indices[i], 2]);
                    res$vertex_region_n_point_z = vappend(res$vertex_region_n_point_z, region_centers_xyz[sorted_region_indices[i], 3]);
                }

            }
        } else {
            stop("Invalid value for 'linkage' parameter.");
        }
    }
    return(data.frame(res));
}


#' @title Append to vector or NULL.
#'
#' @param vec an existing vector or NULL
#'
#' @param val the value to append.
#'
#' @return the value itself if vec is \code{NULL}, the vec with the appended val otherwise.
#'
#' @note There almost certainly is a better way to do this, please open an issue and let me know how if you know it.
#'
#' @keywords internal
vappend <- function(vec, val) {
    if(is.null(vec)) {
        return(val);
    }
    return(c(vec, val));
}


#' @title Find closest regions to coordinate using Euclidean or geodesic distance.
#'
#' @description Finds the closest atlas regions according to the brain surface parcellation 'brainparc' for the given query vertex coordinates. The coordinates must be in the surface space of the surface included in the 'brainparc'.
#'
#' @note This is a wrapper around \code{\link{vertex_closest_regions}}. It first computes the surface vertex closest to the given coordinate (using \code{\link{coord_closest_vertex}}), then runs \code{\link{vertex_closest_regions}} with that vertices coordinates.
#'
#' @inheritParams vertex_closest_regions
#'
#' @inheritParams coord_closest_vertex
#'
#' @examples
#' \dontrun{
#' bp = brainparc_fs(fsbrain::fsaverage.path(), "fsaverage", atlas="aparc");
#' coord_closest_regions(bp, bp$surfaces$white$lh$vertices[1:3,]);
#' }
#'
#'
#' @seealso \code{\link{vertex_closest_regions}} is faster if you have a vertex index for the surface in 'brainparc' instead of a coordinate.
#'
#' @export
coord_closest_regions <- function(brainparc, coordinate, linkage = "single", distance = "euclidean", silent = getOption("brainloc.silent", default = FALSE), num_regions_to_report = 5L) {
    closest_vertex_info = coord_closest_vertex(coordinate, get_surface(brainparc));
    return(vertex_closest_regions(brainparc, closest_vertex_info$both_closest_vertex, closest_vertex_info$both_hemi, linkage = linkage, distance = distance, silent = silent, num_regions_to_report = num_regions_to_report));
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


#' @title Compute centroids of surface atlas parcellation regions.
#'
#' @param bp a brainparc instance
#'
#' @return hemilist of data.frames, the centroid coordinates for the regions of the respective hemispheres.
#'
#' @keywords internal
region_centroids <- function(bp) {
    centr = list();

    regions_lh = vertex_closest_regions(bp, vertices=1L, hemis="lh", linkage = "centroid", distance = "euclidean", num_regions_to_report = NULL);
    centr$lh = regions_lh[c("vertex_region_n_name", "vertex_region_n_point_x", "vertex_region_n_point_y", "vertex_region_n_point_z")];
    colnames(centr$lh) = c("region", "x", "y", "z");

    regions_rh = vertex_closest_regions(bp, vertices=1L, hemis="rh", linkage = "centroid", distance = "euclidean", num_regions_to_report = NULL);
    centr$rh = regions_rh[c("vertex_region_n_name", "vertex_region_n_point_x", "vertex_region_n_point_y", "vertex_region_n_point_z")];
    colnames(centr$rh) = c("region", "x", "y", "z");

    return(centr);
}



#' @title Find closest surface vertex to a point, using Euclidean distance.
#'
#' @param coordinate \code{nx3} numerical matrix or vector of length 3, the query point coordinates.
#'
#' @param surfaces \code{\link{hemilist}} of \code{fs.surface} instances
#'
#' @return a data.frame with columns named 'query_x', 'query_y', 'query_z', 'lh_closest_vertex', 'lh_distance', 'rh_closest_vertex', 'rh_distance', 'both_closest_vertex', 'both_distance', 'both_hemi'.
#'
#' @examples
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
    return(data.frame("query_x"=coordinate[,1], "query_y"=coordinate[,2], "query_z"=coordinate[,3], "lh_closest_vertex"=lh_closest_vertex, "lh_distance"=lh_distance, "rh_closest_vertex"=rh_closest_vertex, "rh_distance"=rh_distance, "both_closest_vertex"=both_closest_vertex, "both_distance"=both_distance, "both_hemi"=both_hemi, stringsAsFactors = FALSE));
}

