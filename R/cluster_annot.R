#' @title Create a hemilist of fs.annot instances from the given cluster overlay.
#'
#' @param clusteroverlay hemilist of integer vectors or a single integer vector of cluster overlay data: a vector that assigns each vertex to an integer class, all vertices belonging to a cluster share the same integer assignment. Non-cluster vertices are assigned the background class, see parameter 'background_code'. One can also pass character strings, which will be interpreted as path to files that should be loaded with \code{freesurferformats::read.fs.morph} to get the clusteroverlay. One can also pass a \code{clusterinfo} instance, the 'overlay' field will be extracted and used in that case.
#'
#' @param background_code scalar integer, the code in the overlayID data that should be interpreted as background or 'not part of any cluster'.
#'
#' @param hemi character string, ignore unless clusteroverlay is a vector. In that case, it must be 'lh' or 'rh'. Used as a prefix when naming the clusters.
#'
#' @return \code{\link{hemilist}} of \code{freesurferformats::fs.annot} instances, each one representing a brain surface parcellation for one hemisphere. Each label in the parcellation (with the exception of the unknown/background region) corresponds to a cluster.
#'
#' @keywords internal
clusteroverlay_to_annot <- function(clusteroverlay, background_code=0L, hemi=NULL) {

    if(is.clusterinfo(clusteroverlay)) { # also accept clusterinfo instances.
        clusteroverlay = clusteroverlay$overlay;
    }

    if(is.list(clusteroverlay)) {
        if(is.character(clusteroverlay$lh)) { # It's a file path, read the file.
            clusteroverlay$lh = freesurferformats::read.fs.morph(clusteroverlay$lh);
        }
        if(is.character(clusteroverlay$rh)) { # It's a file path, read the file.
            clusteroverlay$rh = freesurferformats::read.fs.morph(clusteroverlay$rh);
        }
        annot_lh = clusteroverlay_to_annot(clusteroverlay$lh, background_code=background_code, hemi="lh");
        annot_rh = clusteroverlay_to_annot(clusteroverlay$rh, background_code=background_code, hemi="rh");
        return(list("lh"=annot_lh, "rh"=annot_rh));
    } else {
        if(! (hemi %in% c("lh", "rh"))) {
            stop("If clusteroverlay is a vector, parameter 'hemi' must be one of 'lh' or 'rh'.");
        }
        clusteroverlay = as.integer(clusteroverlay);
        label_vertices_by_region = list();
        region_int_codes = unique(clusteroverlay);

        if(length(region_int_codes) > 1L) { # there are clusters, because there is more than 1 region. If not, its all background.

            current_region_idx = 1L;
            index_of_unknown_region = -1L;
            for(region_code in region_int_codes) {
                if(region_code == background_code) {
                    index_of_unknown_region = current_region_idx;
                } else {
                    label_name = paste(hemi, "cluster", region_code, sep="_");
                    label_vertices_by_region[[label_name]] = which(clusteroverlay == region_code);
                }
                current_region_idx = current_region_idx + 1L;
            }
            if(index_of_unknown_region < 0L) {
                stop(sprintf("No vertex in annot has the background_code '%d', cannot define unknown region.\n", background_code));
            }
        } else {
            index_of_unknown_region = 1L;
            #message(sprintf("Only 1 region code in overlay for hemi '%s': no clusters.\n", hemi));
        }

        num_verts = length(clusteroverlay);
        annot = fsbrain::label.to.annot(label_vertices_by_region, num_vertices_in_surface = num_verts, colortable_df = NULL, index_of_unknown_region = index_of_unknown_region);
        if(! ("fs.annot" %in% class(annot))) {
            class(annot) = c(class(annot), "fs.annot"); # workaround fsbrain bug #36
        }
        return(annot);
    }
}


#' @title Compute number of clusters from cluster annot or clusterinfo instance.
#'
#' @param cluster_annots \code{\link{hemilist}} of cluster annots, see \code{\link{clusteroverlay_to_annot}}.
#'
#' @return named list with keys 'lh', 'rh' and 'total', each holding a scalar integer. The cluster counts.
#'
#' @keywords internal
num_clusters <- function(cluster_annots) {
    if(! freesurferformats::is.fs.annot(cluster_annots$lh)) {
        stop("The 'lh' field of the 'cluster_annots' parameter does not contain a valid fs.annot instance.");
    }
    if(! freesurferformats::is.fs.annot(cluster_annots$rh)) {
        stop("The 'rh' field of the 'cluster_annots' parameter does not contain a valid fs.annot instance.");
    }
    res = list("lh"=0L, "rh"=0L, "total"=0L);
    for (hemi in c("lh", "rh")) {
        for(cluster_name in unique(cluster_annots[[hemi]]$label_names)) {
            if(! (cluster_name %in% c("", "unknown"))) {
                if(hemi == "lh") {
                    res$lh = res$lh + 1L;
                } else {
                    res$rh = res$rh + 1L;
                }
            }
        }
    }
    res$total = res$lh + res$rh;
    return(res);
}


#' @title Given vertex indices defining a cluster, find all atlas regions the cluster overlaps with.
#'
#' @param annot_min a full \code{fs.annot} instance, or a minimal annot, i.e., only the 'label_names' field of the annot. Only a single one, not a \code{\link{hemilist}}.
#'
#' @param cluster_vertices integer vector, the vertices defining the cluster (technically they do not need to form a cluster or be connected). Must not be empty.
#'
#' @return data.frame with columns 'region': the region name, 'num_shared_vertices': the number of cluster vertices which are in the region, and 'percent_shared_vertices': the percent of cluster vertices which are in the region, 'cluster_percent_of_region': how much of the region area is filled by the cluster vertices, in percent. The data.frame rows are ordered descending by 'percent_shared_vertices'.
#'
#' @keywords internal
cluster_overlapping_regions <- function(annot_min, cluster_vertices) {
    if(freesurferformats::is.fs.annot(annot_min)) {
        annot_min = annot_min$label_names;
    }

    if(! is.character(annot_min)) {
        stop("Parameter 'annot_min' must be an fs.annot instance or a vector of character strings.");
    }

    cluster_size = length(cluster_vertices);
    if(cluster_size < 1L) {
        stop("Parameter 'cluster_vertices' must not be empty.");
    }

    # Sanity checks.
    if(any(cluster_vertices > length(annot_min))) {
        stop(sprintf("Annotation has length %d, but maximal vertex index is %d.\n", length(annot_min), max(cluster_vertices)));
    }
    if(any(cluster_vertices < 0L)) {
        stop("All vertex indices in parameter 'cluster_vertices' must be positive integers.");
    }


    overlapping_region_names = unique(annot_min[cluster_vertices]);
    nr = length(overlapping_region_names); # num overlapping regions
    overlapping_region_num_vertex_overlap = rep(0L, nr);
    overlapping_region_percent_overlap = rep(0.0, nr);
    cluster_size_percent_of_region = rep(0.0, nr);

    region_idx = 0L;
    for(region_name in overlapping_region_names) {
        region_idx = region_idx + 1L;
        region_vertex_indices = which(annot_min == region_name);
        overlapping_region_num_vertex_overlap[region_idx] = length(which(annot_min[cluster_vertices] == region_name));
        overlapping_region_percent_overlap[region_idx] = overlapping_region_num_vertex_overlap[region_idx] / cluster_size * 100.0;
        cluster_size_percent_of_region[region_idx] = overlapping_region_num_vertex_overlap[region_idx] / length(region_vertex_indices) * 100.0;
    }
    df = data.frame("region"=overlapping_region_names, "num_shared_vertices"=overlapping_region_num_vertex_overlap, "percent_shared_vertices"=overlapping_region_percent_overlap, "cluster_percent_of_region"=cluster_size_percent_of_region);
    return(df[order(df$percent_shared_vertices),]);
}
