

#' @title Create a hemilist of fs.annot instances from the given cluster overlay.
#'
#' @param clusterdata an integer vector of cluster data: a vector that assigns each vertex to an integer class, all vertices belonging to a cluster share the same integer assignment. Non-cluster vertices are assigned the background class, see parameter 'background_code'.
#'
#' @param background_code scalar integer, the code in the data that should be interpreted as background or 'not part of any cluster'.
#'
#' @param hemi character string, ignore unless clusterdata is a vector. In that case, it must be 'lh' or 'rh'. Used to name clusters, as a prefix.
#'
#' @keywords internal
clusters_to_annot <- function(clusterdata, background_code=0L, hemi=NULL) {

    if(is.list(clusterdata)) {
        annot_lh = clusters_to_annot(clusterdata$lh, background_code=background_code, hemi="lh");
        annot_rh = clusters_to_annot(clusterdata$rh, background_code=background_code, hemi="rh");
        return(list("lh"=annot_lh, "rh"=annot_rh));
    } else {
        if(! (hemi %in% c("lh", "rh"))) {
            stop("If clusterdata is a vector, parameter 'hemi' must be one of 'lh' or 'rh'.");
        }
        clusterdata = as.integer(clusterdata);
        label_vertices_by_region = list();
        region_int_codes = unique(clusterdata);

        current_region_idx = 1L;
        index_of_unknown_region = -1L;
        for(region_code in region_int_codes) {
            if(region_code == background_code) {
                index_of_unknown_region = current_region_idx;
            }
            label_name = paste(hemi, "cluster", region_code, sep="_");
            label_vertices_by_region[[label_name]] = which(clusterdata == region_code);
            current_region_idx = current_region_idx + 1L;
        }
        if(index_of_unknown_region < 0L) {
            stop(sprintf("No vertex in annot has the background_code '%d', cannot define unknown region.\n", background_code));
        }

        num_verts = length(clusterdata);
        return(fsbrain::label.to.annot(label_vertices_by_region, num_vertices_in_surface = num_verts, colortable_df = NULL));
    }
}
