


#' @title Print information on cluster extreme values.
#'
#' @description A cluster extreme value is the most extreme value of a cluster. One can choose whether the minimum, maximum, or absolute maximum is requested. Each cluster has exactly one extreme value (though it may occur at several vertices in rare cases).
#'
#' @param clusterinfo a clusterinfo instance, see the \code{clusterinfo} function to get one.
#'
#' @param type character string, one of 'extreme', 'min' or 'max'. Which cluster value to report. The default of 'extreme' reports the absolutely larger one of the min and the max value and is for the typical use case of \code{t}-value maps.
#'
#' @param silent whether to suppress printing messages to stdout.
#'
#' @param ... passed on to \code{\link{clusteroverlay_to_annot}}.
#'
#' @return a \code{data.frame} with cluster extrema information. The column names should be self-explanatory.
#'
#' @note If the extreme value occurs at several vertices of a cluster, which of these vertices will be reported is undefined.
#'
#' @keywords internal
cluster_extrema <- function(clusterinfo, type = "extreme", silent = getOption("brainloc.silent", default = FALSE), ...) {

    if(! is.clusterinfo(clusterinfo)) {
        stop("Parameter 'clusterinfo' must be a clusterinfo instance.");
    }

    cluster_annots = clusteroverlay_to_annot(clusterinfo$overlay, ...);
    num_clusters_total = num_clusters(cluster_annots)$total;
    if(! silent) {
        cat(sprintf("Computing cluster extrema for %d clusters.\n", num_clusters_total));
    }

    all_cluster_names = rep("?", num_clusters_total);
    all_cluster_hemis = rep("?", num_clusters_total);
    all_cluster_num_vertices = rep(0L, num_clusters_total);
    all_cluster_extreme_value = rep(0.0, num_clusters_total);
    all_cluster_extreme_vertex = rep(0L, num_clusters_total);


    current_cluster_idx = 1L;
    rows_to_remove = c();
    for (hemi in c("lh", "rh")) {
        clusters = get_clusters(clusterinfo, hemi = hemi);
        for(cluster_name in names(clusters)) {
            cluster_vertices = clusters[[cluster_name]];
            cluster_num_vertices = length(cluster_vertices);
            all_cluster_names[current_cluster_idx] = cluster_name;
            all_cluster_hemis[current_cluster_idx] = hemi;
            all_cluster_num_vertices[current_cluster_idx] = cluster_num_vertices;
            cluster_min_statvalue = min(clusterinfo$statmap[[hemi]][cluster_vertices]);
            cluster_max_statvalue = max(clusterinfo$statmap[[hemi]][cluster_vertices]);
            if(is.na(cluster_min_statvalue) | is.na(cluster_max_statvalue)) {
                num_na = sum(is.na(clusterinfo$statmap[[hemi]][cluster_vertices]));
                if(num_na == cluster_num_vertices) {
                    warning(sprintf(" - Cluster #%d %s on hemi %s with %d vertices has %d NA stat data values: min = %s, max = %s. Skipping. Vertex indices are in range %d - %d (%d unique).\n", current_cluster_idx, cluster_name, hemi, cluster_num_vertices, num_na, cluster_min_statvalue, cluster_max_statvalue, min(cluster_vertices), max(cluster_vertices), length(unique(cluster_vertices))));
                    all_cluster_extreme_value[current_cluster_idx] = NA;
                    all_cluster_extreme_vertex[current_cluster_idx] = NA;
                    rows_to_remove = c(rows_to_remove, current_cluster_idx);
                    current_cluster_idx = current_cluster_idx + 1L;
                    next;
                } else {
                    message(sprintf(" - Cluster #%d %s on hemi %s with %d vertices has %d NA stat data values which will be ignored.\n", current_cluster_idx, cluster_name, hemi, cluster_num_vertices, num_na));
                    cluster_min_statvalue = min(clusterinfo$statmap[[hemi]][cluster_vertices], na.rm = TRUE);
                    cluster_max_statvalue = max(clusterinfo$statmap[[hemi]][cluster_vertices], na.rm = TRUE);
                }

            }
            cluster_vertex_with_min_statvalue = cluster_vertices[which.min(clusterinfo$statmap[[hemi]][cluster_vertices])];
            cluster_vertex_with_max_statvalue = cluster_vertices[which.max(clusterinfo$statmap[[hemi]][cluster_vertices])];

            if(type == "extreme") {
                is_pos_greater = abs(cluster_max_statvalue) > abs(cluster_min_statvalue);
                if(is_pos_greater) {
                    cluster_extreme_value = cluster_max_statvalue;
                    cluster_vertex_with_extreme_value = cluster_vertex_with_max_statvalue;
                } else {
                    cluster_extreme_value = cluster_min_statvalue;
                    cluster_vertex_with_extreme_value = cluster_vertex_with_min_statvalue;
                }
            } else if(type == "min") {
                cluster_extreme_value = cluster_min_statvalue;
                cluster_vertex_with_extreme_value = cluster_vertex_with_min_statvalue;
            } else if (type == "max") {
                cluster_extreme_value = cluster_max_statvalue;
                cluster_vertex_with_extreme_value = cluster_vertex_with_max_statvalue;
            } else {
                stop("Invalid 'type' argument. Must be one of 'min', 'max', 'extreme'");
            }
            all_cluster_extreme_value[current_cluster_idx] = cluster_extreme_value;
            all_cluster_extreme_vertex[current_cluster_idx] = cluster_vertex_with_extreme_value;
            if(! silent) {
                cat(sprintf(" - Hemi %s cluster '%s' has size %d vertices and %s stat value %f at vertex %d.\n", hemi, cluster_name, cluster_num_vertices, type, cluster_extreme_value, cluster_vertex_with_extreme_value));
            }

            current_cluster_idx = current_cluster_idx + 1L;
        }

    }

    # Remove rows of clusters with NA values.
    df = data.frame("cluster"=all_cluster_names, "hemi"=all_cluster_hemis, "num_vertices"=all_cluster_num_vertices, "extremum_value"=all_cluster_extreme_value, "extremum_vertex"=all_cluster_extreme_vertex, stringsAsFactors = FALSE);
    if(! is.null(rows_to_remove)) {
        df = df[-rows_to_remove, ];
    }

    return(df);
}


#' @title Compute overlapping atlas regions for clusters.
#'
#' @inheritParams cluster_extrema
#'
#' @export
cluster_region_overlap <- function(clusterinfo, silent = getOption("brainloc.silent", default = FALSE)) {
    if(! is.clusterinfo(clusterinfo)) {
        stop("Parameter 'clusterinfo' must be a clusterinfo instance.");
    }
    if(is.null(clusterinfo$brainparc)) {
        stop("The clusterinfo instance must contain a valid brainparc for this function to be able to work.");
    }

    for(atlas in names(clusterinfo$brainparc$annots)) {
        atlas_annot_min = clusterinfo$brainparc$annots[[atlas]][[hemi]];
        overlap_df = cluster_overlapping_regions(atlas_annot_min, cluster_vertices);
        if(! silent) {
            cat(sprintf("   - Hemi %s cluster '%s' overlaps with %d regions of atlas '%s':\n", hemi, cluster_name, nrow(overlap_df), atlas));
            for(row_idx in seq.int(nrow(overlap_df))) {
                cat(sprintf("     * Region %s: %d of %d cluster vertices in region (%.2f percent). Cluster covers %.2f percent of the region.\n", overlap_df$region[row_idx], overlap_df$num_shared_vertices[row_idx], cluster_num_vertices, overlap_df$percent_shared_vertices[row_idx], overlap_df$cluster_percent_of_region[row_idx]));
            }
        }
    }
    stop("THIS FUNCTION IS WIP");


}


#' @title Compute cluster peaks.
#'
#' @description A peak is a position in a cluster (a vertex) which is assigned the maximum statsmap value over its neighborhood. The neighborhood definition used by this function is the 1-ring neighborhood on the mesh.
#'
#' @inheritParams cluster_extrema
#'
#' @return a \code{data.frame} with cluster peak information. The column names should be self-explanatory.
#'
#' @export
cluster_peaks <- function(clusterinfo, type = "extreme", silent = getOption("brainloc.silent", default = FALSE), ...) {

    if(! is.clusterinfo(clusterinfo)) {
        stop("Parameter 'clusterinfo' must be a clusterinfo instance.");
    }

    cluster_annots = clusteroverlay_to_annot(clusterinfo$overlay, ...);
    num_clusters_total = num_clusters(cluster_annots)$total;

    if(! silent) {
        cat(sprintf("Computing cluster peaks for %d clusters.\n", num_clusters_total));
    }

    all_cluster_names = c();
    all_cluster_hemis = c();
    all_clusters_peak_vertex = c();
    all_clusters_peak_value = c();


    current_cluster_idx = 1L;
    for (hemi in c("lh", "rh")) {
        clusters = get_clusters(clusterinfo, hemi = hemi);
        surface = get_surface(clusterinfo$brainparc)[[hemi]];
        for(cluster_name in names(clusters)) {
            cluster_vertices = clusters[[cluster_name]];
            cluster_num_vertices = length(cluster_vertices);
            cl_peaks = single_cluster_peaks(cluster_vertices, clusterinfo$statmap[[hemi]], surface, type = type);
            num_peaks = length(cl_peaks$vertex); # Could also use cl_peaks$value, the length is identical.
            if(! silent) {
                cat(sprintf("Cluster %s on hemi %s has %d vertices and %d peaks.\n", cluster_name, hemi, length(cluster_vertices), num_peaks));
            }
            if(is.null(all_cluster_names)) {
                all_cluster_names = rep(cluster_name, num_peaks);
                all_cluster_hemis = rep(hemi, num_peaks);
                all_clusters_peak_vertex = cl_peaks$vertex;
                all_clusters_peak_value = cl_peaks$value;
            } else {
                all_cluster_names = c(all_cluster_names, rep(cluster_name, num_peaks));
                all_cluster_hemis = c(all_cluster_hemis, rep(hemi, num_peaks));
                all_clusters_peak_vertex = c(all_clusters_peak_vertex, cl_peaks$vertex);
                all_clusters_peak_value = c(all_clusters_peak_value, cl_peaks$value);
            }
        }
    }
    return(data.frame("cluster"=all_cluster_names, "hemi"=all_cluster_hemis, "peak_vertex"=all_clusters_peak_vertex, "peak_value"=all_clusters_peak_value));
}



#' @title Compute peaks of a single cluster.
#'
#' @inheritParams cluster_extrema
#'
#' @param cluster_vertices integer vector, the vertex indices of the cluster.
#'
#' @param statmap double vector, the full stat map for the surface. Only the values of the cluster_vertices are used.
#'
#' @param surface a single \code{fs.surface} instance. Used to compute the neighborhood of the cluster vertices.
#'
#' @return named list with entries 'vertex' and 'value', they contain an integer vector and a double vector, respectively.
#'
#' @note Called by \code{cluster_peaks}, use that instead.
#'
#' @keywords internal
single_cluster_peaks <- function(cluster_vertices, statmap, surface, type = "extreme") {
    if(! freesurferformats::is.fs.surface(surface)) {
        stop("Parameter 'surface' must be an fs.surface instance.");
    }
    if(nrow(surface$vertices) != length(statmap)) {
        stop(sprintf("Number of surface vertices is %d and does not match length of stat map %d. Lengths must match.\n", nrow(surface$vertices), length(statmap)));
    }
    if(! (is.vector(cluster_vertices) & is.integer(cluster_vertices))) {
        stop("Parameter 'cluster_vertices' must be an integer vector.");
    }
    if(! (is.vector(statmap) & is.numeric(statmap))) {
        stop("Parameter 'statmap' must be an numeric vector.");
    }

    peak_vertices = c();
    peak_values = c();

    adj = Rvcg::vcgVertexNeighbors(fs.surface.to.tmesh3d(surface));
    for(cl_vertex in cluster_vertices) {
        neigh = adj[[cl_vertex]];
        neigh_within_cluster = neigh[neigh %in% cluster_vertices];
        neigh_stat_values = statmap[c(neigh_within_cluster, cl_vertex)];
        cluster_neigh_max_statvalue = max(neigh_stat_values, na.rm = TRUE);
        cluster_neigh_min_statvalue = min(neigh_stat_values, na.rm = TRUE);
        if(type == "extreme") {
            is_pos_greater = abs(cluster_neigh_max_statvalue) > abs(cluster_neigh_min_statvalue);
            if(is_pos_greater) {
                cluster_neigh_extreme_value = cluster_neigh_max_statvalue;
            } else {
                cluster_neigh_extreme_value = cluster_neigh_min_statvalue;
            }
        } else if(type == "min") {
            cluster_neigh_extreme_value = cluster_neigh_min_statvalue;
        } else if (type == "max") {
            cluster_neigh_extreme_value = cluster_neigh_max_statvalue;
        } else {
            stop("Invalid 'type' argument. Must be one of 'min', 'max', 'extreme'");
        }

        if(cluster_neigh_extreme_value == statmap[cl_vertex]) { # The vertex is a peak.
            if(is.null(peak_vertices)) {
                peak_vertices = cl_vertex;
                peak_values = statmap[cl_vertex];
            } else {
                peak_vertices = c(peak_vertices, cl_vertex);
                peak_values = c(peak_values, statmap[cl_vertex]);
            }
        }
    }
    return(list("vertex"=peak_vertices, "value"=peak_values));
}


#' @title Get the clusters from a clusterinfo instance.
#'
#' @param clusterinfo a \code{clusterinfo} instance.
#'
#' @param hemi character string, one of 'lh', 'rh' or 'both'. The hemisphere(s) for which to return the clusters.
#'
#' @return named list, the keys are the cluster names, and the values are integer vectors defining the member vertices.
#'
#' @export
get_clusters <- function(clusterinfo, hemi="both") {
    if(! is.clusterinfo(clusterinfo)) {
        stop("Parameter 'clusterinfo' must be a clusterinfo instance.");
    }
    if(!(hemi %in% c("lh", "rh", "both"))) {
        stop("Parameter 'hemi' must be one of 'lh', 'rh' or 'both'.");
    }
    cluster_annots = clusteroverlay_to_annot(clusterinfo$overlay);

    clusters = list();

    if(hemi == "both") {
        hemi = c("lh", "rh");
    }

    for (hemisphere in hemi) {
        for(cluster_name in unique(cluster_annots[[hemisphere]]$label_names)) {
            if(! (cluster_name %in% c("", "unknown"))) {
                cluster_vertices = which(cluster_annots[[hemisphere]]$label_names == cluster_name);
                clusters[[cluster_name]] = cluster_vertices;
            }
        }
    }
    return(clusters);
}



#' @title Get details on cluster locations.
#'
#' @param clusterinfo a clusterinfo instance, see the \code{clusterinfo} function on how to get one. Must contain a valid \code{brainparc} in field 'brainparc'.
#'
#' @param silent logical, whether to suppress console messages.
#'
#' @param ... passed on to \code{\link{cluster_extrema}}.
#'
#' @return a new version of the input data.frame, with additional columns appended.
#'
#' @export
cluster_location_details <- function(clusterinfo, silent = getOption("brainloc.silent", default = FALSE), ...) {

    if(! is.clusterinfo(clusterinfo)) {
        stop("Parameter 'clusterinfo' must be a clusterinfo instance.");
    }

    extrema = cluster_extrema(clusterinfo, ...);

    if(! silent) {
        cat(sprintf("Computing location details for cluster extrema.\n"));
    }

    if(is.null(clusterinfo$brainparc)) {
        warning("Cannot compute cluster location details, clusterinfo in parameter 'clusterinfo' does not contain a valid 'brainparc'.");
        return(extrema);
    }

    nc = nrow(extrema); # number of clusters

    # results, to be filled.
    all_mni152_R = rep(0.0, nc);
    all_mni152_A = rep(0.0, nc);
    all_mni152_S = rep(0.0, nc);
    all_Tal_R = rep(0.0, nc);
    all_Tal_A = rep(0.0, nc);
    all_Tal_S = rep(0.0, nc);

    for(cluster_idx in seq.int(nc)) {
        hemi = extrema$hemi[cluster_idx];
        surface_hemilist = get_surface(clusterinfo$brainparc);
        surface = surface_hemilist[[hemi]];
        query_vertex = extrema$extremum_vertex[cluster_idx];
        vertex_coords_MNI305 = surface$vertices[query_vertex, ];

        coord_info = coord_MNI305_info(vertex_coords_MNI305);

        all_mni152_R[cluster_idx] = coord_info$mni152[1];
        all_mni152_A[cluster_idx] = coord_info$mni152[2];
        all_mni152_S[cluster_idx] = coord_info$mni152[3];
        all_Tal_R[cluster_idx] = coord_info$talairach[1];
        all_Tal_A[cluster_idx] = coord_info$talairach[2];
        all_Tal_S[cluster_idx] = coord_info$talairach[3];
        if(! silent) {
            cat(sprintf(" - Cluster %s on hemi %s extremum vertex %d has MNI152 coords (%f %f %f) and Talairach coords (%f %f %f).\n", extrema$cluster[cluster_idx], hemi, query_vertex, coord_info$mni152[1], coord_info$mni152[2], coord_info$mni152[3], coord_info$talairach[1], coord_info$talairach[2], coord_info$talairach[3]));
        }
    }

    extrema$mni152_r = all_mni152_R;
    extrema$mni152_a = all_mni152_A;
    extrema$mni152_s = all_mni152_S;
    extrema$talairach_r = all_Tal_R;
    extrema$talairach_a = all_Tal_A;
    extrema$talairach_s = all_Tal_S;
    return(extrema);
}

