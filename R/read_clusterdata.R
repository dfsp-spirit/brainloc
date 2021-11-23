

#' @title Create a hemilist of fs.annot instances from the given cluster overlay.
#'
#' @param clusteroverlay hemilist of integer vectors or a single integer vector of cluster overlay data: a vector that assigns each vertex to an integer class, all vertices belonging to a cluster share the same integer assignment. Non-cluster vertices are assigned the background class, see parameter 'background_code'. One can also pass character strings, which will be interpreted as path to files that should be loaded with \code{freesurferformats::read.fs.morph} to get the clusteroverlay.
#'
#' @param background_code scalar integer, the code in the data that should be interpreted as background or 'not part of any cluster'.
#'
#' @param hemi character string, ignore unless clusteroverlay is a vector. In that case, it must be 'lh' or 'rh'. Used as a prefix when naming the clusters.
#'
#' @return hemilist of \code{freesurferformats::fs.annot} instances, each one representing a brain surface parcellation for one hemisphere. Each label in the parcellation (with the exception of the unknown/background region) corresponds to a cluster.
#'
#' @keywords internal
clusteroverlay_to_annot <- function(clusteroverlay, background_code=1L, hemi=NULL) {

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
    annots = clusteroverlay_to_annot(clusterinfo$overlay);

    num_clusters_total = num_clusters(annots)$total;

    all_cluster_names = rep("?", num_clusters_total);
    all_cluster_hemis = rep("?", num_clusters_total);
    all_cluster_num_vertices = rep(0L, num_clusters_total);
    all_cluster_extreme_value = rep(0.0, num_clusters_total);
    all_cluster_extreme_vertex = rep(0L, num_clusters_total);


    current_cluster_idx = 1L;
    for (hemi in c("lh", "rh")) {
        for(cluster_name in unique(annots[[hemi]]$label_names)) {
            if(! (cluster_name %in% c("", "unknown"))) {
                cluster_vertices = which(annots[[hemi]]$label_names == cluster_name);
                cluster_num_vertices = length(cluster_vertices);
                all_cluster_names[current_cluster_idx] = cluster_name;
                all_cluster_hemis[current_cluster_idx] = hemi;
                all_cluster_num_vertices[current_cluster_idx] = cluster_num_vertices;
                cluster_min_statvalue = min(clusterinfo$statmap[[hemi]][cluster_vertices]);
                cluster_vertex_with_min_statvalue = cluster_vertices[which.min(clusterinfo$statmap[[hemi]][cluster_vertices])];
                cluster_max_statvalue = max(clusterinfo$statmap[[hemi]][cluster_vertices]);
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
                    cat(sprintf("Hemi %s cluster '%s' has size %d vertices and %s stat value %f at vertex %d.\n", hemi, cluster_name, cluster_num_vertices, type, cluster_extreme_value, cluster_vertex_with_extreme_value));
                }
                current_cluster_idx = current_cluster_idx + 1L;
            }
        }
    }
    return(data.frame("cluster"=all_cluster_names, "hemi"=all_cluster_hemis, "num_vertices"=all_cluster_num_vertices, "extremum_value"=all_cluster_extreme_value, "extremum_vertex"=all_cluster_extreme_vertex, stringsAsFactors = FALSE));
}


#' @title Compute number of clusters from cluster annots.
#'
#' @param cluster_annots hemilist of cluster annots, see \code{clusteroverlay_to_annot}.
#'
#' @return named list with keys 'lh', 'rh' and 'total', each holdign a scalar integer. The cluster counts.
#'
#' @keywords internal
num_clusters <- function(cluster_annots) {
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
    res = list("overlay"=list("lh"=lh_overlay, "rh"=rh_overlay), "statmap"=list("lh"=lh_statmap, "rh"=rh_statmap), "metadata"=list("template_subject"=template_subject, "brainparc"=NULL));
    if(dir.exists(file.path(subjects_dir, template_subject))) { # Load surfaces and atlas if possible.
        res$brainparc = brainparc_fs(subjects_dir, template_subject);
    } else {
        message(sprintf("Could not load brainparc for template subject '%s': directory '%s' not found.\n", template_subject, file.path(subjects_dir, template_subject)));
    }

    class(res) = c(class(res), "clusterinfo");
    return(res);
}


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
#' @param clusters a clusterinfo instance, see the \code{clusterinfo} function on how to get one. Must contain a valid \code{brainparc} in field 'brainparc'.
#'
#' @param silent logical, whether to suppress console messages.
#'
#' @return a new version of the input data.frame, with additional columns appended.
get_cluster_location_details <- function(clusters, silent = getOption("brainloc.silent", default = FALSE)) {
    extrema = cluster_extrema(clusters);

    if(is.null(clusters$brainparc)) {
        warning("Cannot compute cluster location details, clusterinfo in parameter 'clusters' does not contain a valid 'brainparc'.");
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
        surface = clusters$brainparc$surfaces$white[[hemi]];
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
            cat(sprintf("Cluster %s extremum vertex %d has MNI152 coords (%f %f %f) and Talairach coords (%f %f %f).\n", extrema$cluster[cluster_idx], query_vertex, coord_info$mni152[1], coord_info$mni152[2], coord_info$mni152[3], coord_info$talairach[1], coord_info$talairach[2], coord_info$talairach[3]));
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


#' @keywords internal
test_clusters_to_annot <- function(sjd = "~/software/freesurfer/subjects", sj="fsaverage") {
    options("brainloc.fs_home"="~/software/freesurfer/");
    lh_an = fsbrain::subject.annot(sjd, sj, hemi="lh", atlas="aparc"); # we abuse an atlas as a cluster overlay file in this example  because it works technically. It does not make any sense, I just did not have a Matlab surfstat output file at hand.
    rh_an = fsbrain::subject.annot(sjd, sj, hemi="rh", atlas="aparc");
    clusteroverlay = list("lh"=strvec2int(lh_an$label_codes), "rh"=strvec2int(rh_an$label_codes));
    thickness = fsbrain::subject.morph.native(sjd, sj, "thickness", hemi="both", split_by_hemi = TRUE);
    tmap = list("lh"=thickness$lh * 2 - 2L, "rh"=thickness$lh * 2 - 2L);  # We abuse a cortical thickness map as a t-value map. Yes, that's really ugly.
    clusters = clusterinfo(clusteroverlay$lh, clusteroverlay$rh, tmap$lh, tmap$rh);
    #cluster_annots = clusteroverlay_to_annot(clusters$overlay);
    #fsbrain::vis.colortable.legend(cluster_annots$lh);
    extrema = cluster_extrema(clusters);
    if(! is.null(clusters$brainparc)) {
        extrema = get_cluster_location_details(clusters);
    }
}
