

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

#' @title Convert vector of character strings to integer vector.
#'
#' @param input vector of character strings
#'
#' @return integer vector
#'
#' @keywords internal
strvec2int <- function(input) { as.integer(as.factor(input)); }


#' @keywords internal
test_clusters_to_annot <- function(sjd = "~/data/tim_only", sj="tim") {
    lh_an = fsbrain::subject.annot(sjd, sj, hemi="lh", atlas="aparc");
    rh_an = fsbrain::subject.annot(sjd, sj, hemi="rh", atlas="aparc");
    clusteroverlay = list("lh"=strvec2int(lh_an$label_codes), "rh"=strvec2int(rh_an$label_codes));
    annots = clusteroverlay_to_annot(clusteroverlay);
    #fsbrain::vis.colortable.legend(annots$lh);
}
