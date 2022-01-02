
#' @title Compute adjacency matrix of surface parcellation (annot) regions.
#'
#' @param annot_min a full \code{fs.annot} instance, or a minimal annot, i.e., only the 'label_names' field of the annot. Only a single one, not a \code{\link{hemilist}}.
#'
#' @param surface \code{fs.surface} instance, the brain mesh (for a single hemisphere).
#'
#' @return integer matrix of regions, expressing whether they are direct neighbors (value \code{1L}) or not (value \code{0L}).
#'
#' @note Regions from different hemispheres cannot be direct neighbors.
#'
#' @note This works on a single hemisphere. We should create a wrapper that works on brainparc instances.
#'
#' @keywords internal
annot_neighbors <- function(annot_min, surface) {
    if(freesurferformats::is.fs.annot(annot_min)) {
        annot_min = get_minannot(annot_min);
    }

    if(! is.character(annot_min)) {
        stop("Parameter 'annot_min' must be an fs.annot instance or a vector of character strings.");
    }

    region_names = unique(annot_min);
    nr = length(region_names); # num regions

    region_adj = matrix(data = rep(0L, (nr*nr)), ncol = nr);
    adj = Rvcg::vcgVertexNeighbors(fs.surface.to.tmesh3d(surface));

}
