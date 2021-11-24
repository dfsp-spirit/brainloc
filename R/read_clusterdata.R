

#' @title Create a hemilist of fs.annot instances from the given cluster overlay.
#'
#' @param clusteroverlay hemilist of integer vectors or a single integer vector of cluster overlay data: a vector that assigns each vertex to an integer class, all vertices belonging to a cluster share the same integer assignment. Non-cluster vertices are assigned the background class, see parameter 'background_code'. One can also pass character strings, which will be interpreted as path to files that should be loaded with \code{freesurferformats::read.fs.morph} to get the clusteroverlay. One can also pass a \code{clusterinfo} instance, the 'overlay' field will be extracted and used in that case.
#'
#' @param background_code scalar integer, the code in the data that should be interpreted as background or 'not part of any cluster'.
#'
#' @param hemi character string, ignore unless clusteroverlay is a vector. In that case, it must be 'lh' or 'rh'. Used as a prefix when naming the clusters.
#'
#' @return hemilist of \code{freesurferformats::fs.annot} instances, each one representing a brain surface parcellation for one hemisphere. Each label in the parcellation (with the exception of the unknown/background region) corresponds to a cluster.
#'
#' @keywords internal
clusteroverlay_to_annot <- function(clusteroverlay, background_code=1L, hemi=NULL) {

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

        num_verts = length(clusteroverlay);
        annot = fsbrain::label.to.annot(label_vertices_by_region, num_vertices_in_surface = num_verts, colortable_df = NULL, index_of_unknown_region = index_of_unknown_region);
        if(! ("fs.annot" %in% class(annot))) {
            class(annot) = c(class(annot), "fs.annot"); # workaround fsbrain bug #36
        }
        return(annot);
    }
}


#' @title Print information on cluster extreme values.
#'
#' @param clusterinfo a clusterinfo instance, see the \code{clusterinfo} function to get one.
#'
#' @param type character string, one of 'extreme', 'min' or 'max'. Which cluster value to report. The default of 'extreme' reports the absolutely larger one of the min and the max value and is for the typical use case of \code{t}-value maps.
#'
#' @param silent whether to suppress printing messages to stdout.
#'
#' @return a \code{data.frame} with cluster information. The column names should be self-explanatory.
#'
#' @keywords internal
cluster_extrema <- function(clusterinfo, type = "extreme", silent = getOption("brainloc.silent", default = FALSE)) {

    if(! is.clusterinfo(clusterinfo)) {
        stop("Parameter 'clusterinfo' must be a clusterinfo instance.");
    }

    if(! silent) {
        cat(sprintf("Computing cluster extrema.\n"));
    }

    cluster_annots = clusteroverlay_to_annot(clusterinfo$overlay);

    num_clusters_total = num_clusters(cluster_annots)$total;

    all_cluster_names = rep("?", num_clusters_total);
    all_cluster_hemis = rep("?", num_clusters_total);
    all_cluster_num_vertices = rep(0L, num_clusters_total);
    all_cluster_extreme_value = rep(0.0, num_clusters_total);
    all_cluster_extreme_vertex = rep(0L, num_clusters_total);


    current_cluster_idx = 1L;
    rows_to_remove = c();
    for (hemi in c("lh", "rh")) {
        for(cluster_name in unique(cluster_annots[[hemi]]$label_names)) {
            if(! (cluster_name %in% c("", "unknown"))) {
                cluster_vertices = which(cluster_annots[[hemi]]$label_names == cluster_name);
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


                if(! is.null(clusterinfo$brainparc)) {
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
                }


                current_cluster_idx = current_cluster_idx + 1L;
            }
        }
    }

    # Remove rows of clusters with NA values.
    df = data.frame("cluster"=all_cluster_names, "hemi"=all_cluster_hemis, "num_vertices"=all_cluster_num_vertices, "extremum_value"=all_cluster_extreme_value, "extremum_vertex"=all_cluster_extreme_vertex, stringsAsFactors = FALSE);
    if(! is.null(rows_to_remove)) {
        df = df[-rows_to_remove, ];
    }

    return(df);
}


#' @title Compute number of clusters from cluster annot or clusterinfo instance.
#'
#' @param cluster_annots hemilist of cluster annots, see \code{\link{clusteroverlay_to_annot}}.
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


#' @title Create clusterinfo data structure from files or pre-loaded data.
#'
#' @param lh_overlay integer vector, the cluster overlay data: one integer per vertex of the left hemisphere. This assigns each vertex to a cluster, and all vertices of one cluster have the same number. Background is typically 1L. If a character string, the parameter will be interpreted as file path and loaded with \code{freesurferformats::read.fs.morph}.
#'
#' @param rh_overlay integer vector, just like \code{lh_overlay}, but for the right hemisphere.
#'
#' @param lh_statmap double vector, the stats map. Typically a t-value map. Must have one value per vertex. If a character string, the parameter will be interpreted as file path and loaded with \code{freesurferformats::read.fs.morph}.
#'
#' @param rh_statmap double vector, just like \code{lh_statmap}, but for the right hemisphere.
#'
#' @param template_subject character string, the template subject name. Typically 'fsaverage' or 'fsaverage6'. Must be some MNI305 space subject.
#'
#' @param subjects_dir character string, file system path to a directory containing the recon-all data for the template_subject. Used to load surfaces and annotations to identify cluster coordinates and atlas regions.
#'
#' @return named list with entries 'overlay', 'statmap', and 'metadata': a clusterinfo data structure. Each of the 'overlay' and 'statmap' keys holds a hemilist of numerical vectors.
#'
#' @export
clusterinfo <- function(lh_overlay, rh_overlay, lh_statmap, rh_statmap, template_subject="fsaverage", subjects_dir=file.path(getOption("brainloc.fs_home", default = Sys.getenv("FREESURFER_HOME")), 'subjects')) {
    if(is.character(lh_overlay)) {
        lh_overlay = freesurferformats::read.fs.morph(lh_overlay);
    }
    if(is.character(rh_overlay)) {
        rh_overlay = freesurferformats::read.fs.morph(rh_overlay);
    }
    if(is.character(lh_statmap)) {
        lh_statmap = freesurferformats::read.fs.morph(lh_statmap);
    }
    if(is.character(rh_statmap)) {
        rh_statmap = freesurferformats::read.fs.morph(rh_statmap);
    }

    # Some sanity checks.
    if(length(lh_overlay) != length(lh_statmap)) {
        stop(sprintf("Data invalid: left hemisphere overlay and statmap must have identical lengths, but have %d versus %d.\n", length(lh_overlay), length(lh_statmap)));
    }
    if(length(rh_overlay) != length(rh_statmap)) {
        stop(sprintf("Data invalid: right hemisphere overlay and statmap must have identical lengths, but have %d versus %d.\n", length(rh_overlay), length(rh_statmap)));
    }

    res = list("overlay"=list("lh"=lh_overlay, "rh"=rh_overlay), "statmap"=list("lh"=lh_statmap, "rh"=rh_statmap), "metadata"=list("template_subject"=template_subject, "brainparc"=NULL));
    if(dir.exists(file.path(subjects_dir, template_subject))) { # Load surfaces and atlas if possible.
        res$brainparc = brainparc_fs(subjects_dir, template_subject);
        # This allows us to do more sanity checks.
        # Some sanity checks.
        lh_surface_numverts = nrow(res$brainparc$surfaces[[1]]$lh$vertices);
        if(length(lh_overlay) != lh_surface_numverts) {
            stop(sprintf("Data invalid: left hemisphere overlay size must match surface vertex count, but they differ: %d versus %d.\n", length(lh_overlay), lh_surface_numverts));
        }
        rh_surface_numverts = nrow(res$brainparc$surfaces[[1]]$rh$vertices);
        if(length(rh_overlay) != rh_surface_numverts) {
            stop(sprintf("Data invalid: right hemisphere overlay size must match surface vertex count, but they differ: %d versus %d.\n", length(rh_overlay), rh_surface_numverts));
        }
    } else {
        message(sprintf("Could not load brainparc for template subject '%s': directory '%s' not found.\n", template_subject, file.path(subjects_dir, template_subject)));
    }

    class(res) = c(class(res), "clusterinfo");
    return(res);
}


#' @title Get the clusters from a clusterinfo instance.
#'
#' @param clusterinfo a \code{clusterinfo} instance.
#'
#' @return named list, the keys are the cluster names, and the values are integer vectors defining the member vertices.
#'
#' @keywords internal
get_clusters <- function(clusterinfo) {
    if(! is.clusterinfo(clusterinfo)) {
        stop("Parameter 'clusterinfo' must be a clusterinfo instance.");
    }
    cluster_annots = clusteroverlay_to_annot(clusterinfo$overlay);

    clusters = list();

    for (hemi in c("lh", "rh")) {
        for(cluster_name in unique(cluster_annots[[hemi]]$label_names)) {
            if(! (cluster_name %in% c("", "unknown"))) {
                cluster_vertices = which(cluster_annots[[hemi]]$label_names == cluster_name);
                clusters[[cluster_name]] = cluster_vertices;
            }
        }
    }
    return(clusters);
}


#' @title Check whether x is a clusterinfo instance.
#'
#' @param x any R object
#'
#' @return logical, whether x is a clusterinfo instance.
#'
#' @export
is.clusterinfo <- function(x) inherits(x, "clusterinfo")


#' @title Convert vector of character strings to integer vector.
#'
#' @param input vector of character strings
#'
#' @return integer vector
#'
#' @keywords internal
strvec2int <- function(input) { as.integer(as.factor(input)); }


#' @title Get details on cluster location, required brainparc.
#'
#' @param clusterinfo a clusterinfo instance, see the \code{clusterinfo} function on how to get one. Must contain a valid \code{brainparc} in field 'brainparc'.
#'
#' @param silent logical, whether to suppress console messages.
#'
#' @return a new version of the input data.frame, with additional columns appended.
get_cluster_location_details <- function(clusterinfo, silent = getOption("brainloc.silent", default = FALSE)) {

    if(! is.clusterinfo(clusterinfo)) {
        stop("Parameter 'clusterinfo' must be a clusterinfo instance.");
    }

    extrema = cluster_extrema(clusterinfo);

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


#' @title Internal test function, will be gone soon and converted into separate unit tests.
#'
#' @note This function is NOT part of the API, using it in client code is a programming error.
#'
#' @keywords internal
test_clusters_to_annot <- function(sjd = "~/software/freesurfer/subjects", sj="fsaverage") {
    options("brainloc.fs_home"="~/software/freesurfer/");
    lh_an = fsbrain::subject.annot(sjd, sj, hemi="lh", atlas="aparc"); # we abuse an atlas as a cluster overlay file in this example  because it works technically. It does not make any sense, I just did not have a Matlab surfstat output file at hand.
    rh_an = fsbrain::subject.annot(sjd, sj, hemi="rh", atlas="aparc");
    clusteroverlay = list("lh" = strvec2int(lh_an$label_codes), "rh" = strvec2int(rh_an$label_codes));
    thickness = fsbrain::subject.morph.native(sjd, sj, "thickness", hemi="both", split_by_hemi = TRUE);
    tmap = list("lh"=thickness$lh * 2 - 2L, "rh"=thickness$lh * 2 - 2L);  # We abuse a cortical thickness map as a t-value map. Yes, that's really ugly.
    clinfo = clusterinfo(clusteroverlay$lh, clusteroverlay$rh, tmap$lh, tmap$rh);
    #cluster_annots = clusteroverlay_to_annot(clinfo$overlay);
    #fsbrain::vis.colortable.legend(cluster_annots$lh);
    extrema_details = get_cluster_location_details(clinfo);
}


#' @title Given vertex indices defining a cluster, find all atlas regions the cluster overlaps with.
#'
#' @param annot_min a full \code{fs.annot} instance, or a minimal annot, i.e., only the 'label_names' field of the annot. Only a single one, not a hemilist.
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
