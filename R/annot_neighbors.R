
#' @title Compute adjacency matrix of surface parcellation (annot) regions.
#'
#' @description Given a brain surface and a parcellation into regions, find out which regions are adjacent to each other. A region \code{i} is adjacent to another region \code{j} if an edge connects any vertex of
#'\code{i} with any vertex of \code{i}.
#'
#' @param annot_min a full \code{fs.annot} instance, or a minimal annot, i.e., only the 'label_names' field of the annot. Only a single one, not a \code{\link{hemilist}}.
#'
#' @param surface \code{fs.surface} instance, the brain mesh (for a single hemisphere).
#'
#' @param empty_rename optional character string, the name to use for unnamed regions (which have an empty string as region name).
#'
#' @return named integer matrix of regions, expressing whether they are direct neighbors (value \code{1L}) or not (value \code{0L}).
#'
#' @note Regions from different hemispheres cannot be direct neighbors.
#'
#' @note This works on a single hemisphere. We should create a wrapper that works on brainparc instances.
#'
#' @keywords internal
annot_neighbors <- function(annot_min, surface, empty_rename="_") {
    if(freesurferformats::is.fs.annot(annot_min)) {
        annot_min = get_minannot(annot_min);
    }

    if(! is.character(annot_min)) {
        stop("Parameter 'annot_min' must be an fs.annot instance or a vector of character strings.");
    }

    if( ! freesurferformats::is.fs.surface(surface)) {
        stop("Parameter 'surface' must be an fs.surface instance.");
    }

    region_names = sort(unique(annot_min));

    # Fix empty region name, which would lead to errors below with the matrix indexing by names.
    region_names[nchar(region_names) == 0] = empty_rename;
    annot_min[nchar(annot_min) == 0] = empty_rename;


    region_mtx_indices = as.integer(as.factor(sort(unique(region_names)))); # this is an arbitrary, local index, based on alphabetical order of the region names. It is used to define the matrix index of a region.
    nr = length(region_names); # num regions
    nv = nrow(surface$vertices); # num vertices

    if(nv != length(annot_min)) {
        stop(sprintf("Surface has %d vertices, but annot_min is for %d vertices. Lengths must match.\n", nv, length(annot_min)));
    }

    region_adj = diag(nr);
    colnames(region_adj) = region_names;
    rownames(region_adj) = region_names;
    adj = Rvcg::vcgVertexNeighbors(fs.surface.to.tmesh3d(surface));

    for(vidx in seq.int(nv)) {
        src_region = annot_min[vidx];
        neigh_verts = adj[[vidx]];

        for(neigh_reg in unique(annot_min[neigh_verts])) {
            #cat(sprintf("Setting adjacency between regions '%s' and '%s'.\n", src_region, neigh_reg));
            region_adj[src_region, neigh_reg] = 1L;
        }
    }
    return(region_adj);
}


#' @title Compute adjacency of surface parcellation (annot) regions.
#'
#' @description Given a brain surface and a parcellation into regions, find out which regions are adjacent to each other. A region \code{i} is adjacent to another region \code{j} if an edge connects any vertex of
#'\code{i} with any vertex of \code{i}. Works with brainparcellations.
#'
#' @param bp a brainparc instance
#'
#' @return hemilist of named integer matrix of regions, expressing whether they are direct neighbors (value \code{1L}) or not (value \code{0L}).
#'
#' @export
brainparc_neighbors <- function(bp) {
    if(! is.brainparc(bp)) {
        stop("Parameter 'bp' must be a brainparc instance.");
    }
    bp_neigh = list();
    for(atlas in names(bp$annots)) {
        bp_neigh[[atlas]] = list("lh"=annot_neighbors(bp$annots[[atlas]]$lh, get_surface(bp)$lh), "rh"=annot_neighbors(bp$annots[[atlas]]$rh, get_surface(bp)$rh));
    }
    return(bp_neigh);
}

