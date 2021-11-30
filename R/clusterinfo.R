# Functions for parsing files/pvd and creating clusterinfo instances from them.

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
#' @return named list with entries 'overlay', 'statmap', and 'metadata': a clusterinfo data structure. Each of the 'overlay' and 'statmap' keys holds a \code{\link{hemilist}} of numerical vectors.
#'
#' @seealso \code{\link{clusterinfo_from_thresholded_overlay}} If you do have a thresholded t-map instead of one t-map and one cluster overlay map.
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

    res = list("overlay"=list("lh"=as.integer(lh_overlay), "rh"=as.integer(rh_overlay)), "statmap"=list("lh"=lh_statmap, "rh"=rh_statmap), "metadata"=list("template_subject"=template_subject));
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


#' @title Create clusterinfo instance from thresholded per-vertex data overlays.
#'
#' @inheritParams clusterinfo
#'
#' @param lh_threshmap double vector, the stats map. Typically a thresholded t-value map. Must have one value per vertex. The value assigned to vertices that have been removed by the thresholding can be set with parameter 'value_thresholded'. If a character string, the parameter will be interpreted as file path and loaded with \code{freesurferformats::read.fs.morph}.
#'
#' @param rh_threshmap like \code{lh_threshmap}, but for the right hemisphere.
#'
#' @param value_thresholded scalar double, the data value of thresholded vertices in the \code{lh_threshmap} and \code{rh_threshmap}.
#'
#' @note This can only work if there are no clusters which touch each other. I think this can never happen, but we should double-check.
#'
#' @return a \code{clusterinfo} instance.
#'
#' @seealso \code{\link{clusterinfo}} If you do have a separate t-map and a cluster overlay map instead of only a thresholded t-map.
#'
#' @export
clusterinfo_from_thresholded_overlay <- function(lh_threshmap, rh_threshmap, value_thresholded = 0.0, template_subject="fsaverage", subjects_dir=file.path(getOption("brainloc.fs_home", default = Sys.getenv("FREESURFER_HOME")), 'subjects')) {
    surfaces = subject.surface(subjects_dir, template_subject, surface = "white");
    lh_overlay = clusteroverlay_from_threshmap(lh_threshmap, value_thresholded = value_thresholded, surface = surfaces$lh);
    rh_overlay = clusteroverlay_from_threshmap(rh_threshmap, value_thresholded = value_thresholded, surface = surfaces$rh);
    return(clusterinfo(lh_overlay, rh_overlay, lh_threshmap, rh_threshmap, template_subject = template_subject, subjects_dir = subjects_dir));
}


#' @title Compute overlayID map from thresholded statmap (clustermap).
#'
#' @description This uses breadth-first search on the graph of the mesh to identify all connected cluster values and identify the separate clusters. Different clusters must not touch for this to work (if they do, there is no reliable way to obtain the clusters from the thresholded map). Clusters can have holes.
#'
#' @inheritParams clusterinfo_from_thresholded_overlay
#'
#' @param theshmap double vector, the stats map. Typically a thresholded t-value map (cluster map). Must have one value per vertex. The value assigned to vertices that have been removed by the thresholding can be set with parameter 'value_thresholded'.
#'
#' @param surface a single \code{fs.surface} instance, used for neighborhood computation.
#'
#' @return integer vector, the clusterID overlay.
#'
#' @note This visits every node and edge in the mesh/graph exactly once.
#'
#' @keywords internal
clusteroverlay_from_threshmap <- function(threshmap, value_thresholded, surface) {
    if(! freesurferformats::is.fs.surface(surface)) {
        stop("Parameter 'surface' must be an fs.surface instance.");
    }

    adj = Rvcg::vcgVertexNeighbors(fs.surface.to.tmesh3d(surface));
    num_vertices = nrow(surface$vertices);
    if(length(threshmap) != num_vertices) {
        stop(sprintf("Thresholded map has %d values, but surface has %d vertices. Counts must match.\n", length(threshmap), num_vertices));
    }
    vertex_visited = rep(FALSE, num_vertices);
    overlay_background_value = 0L;
    overlay_unset_value = -1L; # not known yet.
    overlay = rep(overlay_unset_value, num_vertices);

    vertex_visited[threshmap == value_thresholded] = TRUE;
    overlay[threshmap == value_thresholded] = overlay_background_value;
    current_cluster_label_int = overlay_background_value;
    for(start_vertex in seq.int(num_vertices)) {
        if(vertex_visited[start_vertex]) {    # Start BFS at every foreground/cluster vertex and mark all neighbors that can be reached without crossing background vertices.
            next;
        }
        current_cluster_label_int = current_cluster_label_int + 1L;   # Start a new cluster: when we get here this vertex was not reachable from the previous one.
        q = dequer::queue();
        dequer::pushback(q, start_vertex);
        while(length(q) > 0L) {
            v = dequer::pop(q);
            vertex_visited[v] = TRUE;
            if(threshmap[v] != value_thresholded) {
                overlay[v] = current_cluster_label_int;
            }
            for(v_neighbor in adj[[v]]) {
                if(! vertex_visited[v_neighbor]) {
                    vertex_visited[v_neighbor] = TRUE;
                    dequer::pushback(q, v_neighbor);
                }
            }
        }
    }
    return(overlay);
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


