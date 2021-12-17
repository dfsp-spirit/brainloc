

#' @title Create a brain parcellation from a subject in a FreeSurfer recon-all output directory
#'
#' @inheritParams subject.surface
#'
#' @inheritParams subject.annot
#'
#' @param atlas vector of character strings, the brain atlases to use. For this FreeSurfer version, something like "aparc" for the Desikan-Killiani atlas, "aparc.a2009s" for Destrieux, etc.
#'
#' @return a \code{brainparc} instance
#'
#' @export
brainparc_fs <- function(subjects_dir, subject_id, surface="white", atlas=c("aparc")) {
    surfaces_hemilist = subject.surface(subjects_dir, subject_id, hemi="both", surface = surface);
    surfaces = list();
    surfaces[[surface]] = surfaces_hemilist;
    annots = list();
    for(atlas_name in atlas) {
        annots[[atlas_name]] = subject_minannot_fs(subjects_dir, subject_id, atlas=atlas_name);
    }
    return(brainparc(surfaces, annots));
}


#' @title Create a brain parcellation from custom data.
#'
#' @param surfaces a named list with a single entry, which must contain a hemilist of fs.surface instances. The name can be anything describing the surface. Typically it is something like 'white', 'pial', or 'inflated' for FreeSurfer surfaces.
#'
#' @param annots a named list with at least one entry, all of which must contain a hemilist of minannot (or \code{fs.annot}) instances. The name can be anything describing the atlas from which the parcellation originates. Typically it is something like 'aparc', 'aparc.a2009s', or 'Desikan' for FreeSurfer parcellations.
#'
#' @note A \code{minannot} for a surface is just a vector with length equal to the number of surface vertices, that holds one area label per vertex. An area label is an arbitrary character string, but typically a brain region name.
#'
#' @return a \code{brainparc} instance
#'
#' @export
brainparc <- function(surfaces, annots) {
    if(! is.list(surfaces)) {
        stop("Parameter 'surfaces' must be a list of containing a single hemilist of surfaces, but given one is not a list.");
    }
    if(length(surfaces) != 1L) {
        stop("Parameter 'surfaces' must be a list of containing a single hemilist of surfaces, but given list does not have length 1.");
    }
    if(! is.list(surfaces[[1]])) {
        stop("Parameter 'surfaces' must be a list of containing a single hemilist of surfaces, but first entry of outer list is not a (hemi)list.");
    }
    if(! is.list(annots)) {
        stop("Parameter 'annots' must be a list of containing a single hemilist of minannots, but given one is not a list.");
    }
    if(length(annots) < 1L) {
        stop("Parameter 'annots' must be a list of containing min annots for at least 1 atlas, but list is empty.");
    }
    for(atlas_name in names(annots)) {
        if(! is.list(annots[[atlas_name]])) {
            stop(sprintf("Parameter 'annots' must be a named list containing one or more atlas parcellation hemilists, but entry '%s' of outer list is not a (hemi)list.\n", atlas_name));
        }
        annots[[atlas_name]] = get_minannot(annots[[atlas_name]]); # Convert possible fs.annot to minannot, fails if cannot be done.
    }
    ret = list("surfaces"=surfaces, "annots"=annots);
    class(ret) = c(class(ret), "brainparc");
    return(ret);
}


#' @title Extract the surface hemilist from a brainparc.
#'
#' @param brainparc a \code{brainparc} instance, see \code{\link{brainparc_fs}} to get one.
#'
#' @return \code{\link{hemilist}} of \code{fs.surface} instances.
#'
#' @family brainparc accessors
#'
#' @export
get_surface <- function(brainparc) {
    if(! is.brainparc(brainparc)) {
        stop("Parameter 'brainparc' must be a brainparc instance.");
    }
    return(brainparc$surfaces[[1]]);
}


#' @title Extract the vertex coordinates of the brainparc surface.
#'
#' @param brainparc a \code{brainparc} instance, see \code{\link{brainparc_fs}} to get one.
#'
#' @param vertices integer vector, the query vertex indices.
#'
#' @param hemis vector of character string, the hemispheres for the query vertices. Length must match, length of 1 will be recycled to all vertices.
#'
#' @return numerical nx3 matrix, the coordindates for the x query vertices.
#'
#' @family brainparc accessors
#'
#' @export
get_surface_coords <- function(brainparc, vertices, hemis) {
    nv = length(vertices);
    if(nv < 1L) {
        stop("Parameter 'vertices' must not be empty.");
    }
    if(nv !=  length(hemis)) {
        if(length(hemis) == 1L) {
            hemis = rep(hemis, length(vertices));
        } else {
            stop("Parameter 'hemis' must have length equal to 'vertices' parameter, or exactly length 1.");
        }
    }
    if(any(!(hemis %in% c("lh", "rh")) )) {
        stop("All entries of parameter 'hemis' must be one of 'lh' or 'rh'.");
    }
    surfaces = get_surface(brainparc);
    coords = matrix(rep(0.0, (nv*3L)), nrow = nv);
    for(vidx in seq_along(vertices)) {
        hemi = hemis[[vidx]];
        surf_vertex_idx = vertices[vidx];
        coords[vidx, ] = surfaces[[hemi]]$vertices[surf_vertex_idx, ];
    }
    return(coords);
}


#' @title Return simplified annotation: the region names only.
#'
#' @inheritParams subject.annot
#'
#' @keywords internal
subject_minannot_fs <- function(subjects_dir, subject_id, atlas) {
    lh_annot = subject.annot(subjects_dir, subject_id, hemi="lh", atlas);
    rh_annot = subject.annot(subjects_dir, subject_id, hemi="rh", atlas);
    return(list("lh"=get_minannot(lh_annot), "rh"=get_minannot(rh_annot)));
}


#' @title Extract minannot from fs.annot instance.
#'
#' @param fs.annot a \code{freesurferformats::fs.annot} instance, or a \code{minannot} instance, which will be returned as-is. Alternatively hemilist of the latter.
#'
#' @return a \code{minannot} instance, i.e., only the 'label_names' field of the \code{fs.annot} instance.
#'
#' @keywords internal
get_minannot <- function(fs.annot) {
    if(is.hemilist(fs.annot)) {
        res = list();
        if("lh" %in% names(fs.annot)) {
            res$lh = get_minannot(fs.annot$lh);
        }
        if("rh" %in% names(fs.annot)) {
            res$rh = get_minannot(fs.annot$rh);
        }
        return(res);
    }
    if(is.minannot(fs.annot)) {
        return(fs.annot);
    }
    if(! freesurferformats::is.fs.annot(fs.annot)) {
        stop("Parameter 'fs.annot' must be an fs.annot or minannot instance.");
    }
    minannot = fs.annot$label_names;
    class(minannot) = c(class(minannot), "minannot");
    return(minannot);
}


#' @title Check whether object is a minannot instance.
#'
#' @param x any `R` object
#'
#' @return logical, \code{TRUE} if parameter \code{x} is a minannot instance (that is, has "minannot" in its classes) and \code{FALSE} otherwise.
#'
#' @keywords internal
is.minannot <- function(x) inherits(x, "minannot")



#'@title Load an annotation for a subject.
#'
#' @description Load a brain surface annotation, i.e., a cortical parcellation based on an atlas, for a subject.
#'
#' @param subjects_dir string. The FreeSurfer SUBJECTS_DIR, i.e., a directory containing the data for all your subjects, each in a subdir named after the subject identifier.
#'
#' @param subject_id string. The subject identifier
#'
#' @param hemi string, one of 'lh' or 'rh'. The hemisphere name. Used to construct the names of the annotation and morphometry data files to be loaded.
#'
#' @param atlas string. The atlas name. E.g., "aparc", "aparc.2009s", or "aparc.DKTatlas". Used to construct the name of the annotation file to be loaded.
#'
#' @return the annotation, as returned by \code{\link[freesurferformats]{read.fs.annot}}. It is a named list, enties are: "vertices" vector of n vertex indices, starting with 0. "label_codes": vector of n integers, each entry is a color code, i.e., a value from the 5th column in the table structure included in the "colortable" entry (see below). "label_names": the n brain structure names for the vertices, already retrieved from the colortable using the code. "hex_colors_rgb": Vector of hex color for each vertex.
#'      The "colortable" is another named list with 3 entries: "num_entries": int, number of brain structures. "struct_names": vector of strings, the brain structure names. "table": numeric matrix with num_entries rows and 5 colums. The 5 columns are: 1 = color red channel, 2=color blue channel, 3=color green channel, 4=color alpha channel, 5=unique color code. "colortable_df": The same information as a dataframe. Contains the extra columns "hex_color_string_rgb" and "hex_color_string_rgba" that hold the color as an RGB(A) hex string, like "#rrggbbaa".
#'
#' @family atlas functions
#'
#' @keywords internal
subject.annot <- function(subjects_dir, subject_id, hemi, atlas) {

    if(!(hemi %in% c("lh", "rh", "both"))) {
        stop(sprintf("Parameter 'hemi' must be one of 'lh', 'rh' or 'both' but is '%s'.\n", hemi));
    }

    if(hemi == "both") {
        lh_annot_file = file.path(subjects_dir, subject_id, "label", sprintf("%s.%s.annot", "lh", atlas));
        if(!file.exists(lh_annot_file)) {
            stop(sprintf("Annotation lh file '%s' for subject '%s' atlas '%s' hemi '%s' cannot be accessed.\n", lh_annot_file, subject_id, atlas, "lh"));
        }
        lh_annot = freesurferformats::read.fs.annot(lh_annot_file);

        rh_annot_file = file.path(subjects_dir, subject_id, "label", sprintf("%s.%s.annot", "rh", atlas));
        if(!file.exists(rh_annot_file)) {
            stop(sprintf("Annotation rh file '%s' for subject '%s' atlas '%s' hemi '%s' cannot be accessed.\n", rh_annot_file, subject_id, atlas, "rh"));
        }
        rh_annot = freesurferformats::read.fs.annot(rh_annot_file);

        merged_annot = merge.hemi.annots(lh_annot, rh_annot);
        return(merged_annot);
    }
    else {
        annot_file = file.path(subjects_dir, subject_id, "label", sprintf("%s.%s.annot", hemi, atlas));
        if(!file.exists(annot_file)) {
            stop(sprintf("Annotation file '%s' for subject '%s' atlas '%s' hemi '%s' cannot be accessed.\n", annot_file, subject_id, atlas, hemi));
        }
        return(freesurferformats::read.fs.annot(annot_file));
    }
}


#'@title Merge the annotations from two hemispheres into one annot.
#
#' @param lh_annot, annot. An annotation, as returned by freesurferformats::read.fs.annot().
#'
#' @param rh_annot, annot. An annotation, as returned by freesurferformats::read.fs.annot().
#'
#' @return annot, the merged annotation.
#'
#' @keywords internal
merge.hemi.annots <- function(lh_annot, rh_annot) {
    merged_annot = list();
    merged_annot$colortable = lh_annot$colortable;        # randomly use the lh one, they must be identical for lh nad rh anyways
    merged_annot$colortable_df = lh_annot$colortable_df;  # same

    merged_annot$vertices = c(lh_annot$vertices, rh_annot$vertices);
    merged_annot$label_codes = c(lh_annot$label_codes, rh_annot$label_codes);
    merged_annot$label_names = c(lh_annot$label_names, rh_annot$label_names);
    merged_annot$hex_colors_rgb = c(lh_annot$hex_colors_rgb, rh_annot$hex_colors_rgb);
    return(merged_annot);
}


#' @title Load a surface for a subject.
#'
#' @description Load a brain surface mesh for a subject.
#'
#' @param subjects_dir string. The FreeSurfer SUBJECTS_DIR, i.e., a directory containing the data for all your subjects, each in a subdir named after the subject identifier.
#'
#' @param subject_id string. The subject identifier
#'
#' @param surface string. The surface name. E.g., "white", or "pial". Used to construct the name of the surface file to be loaded.
#'
#' @param hemi string, one of 'lh', 'rh', or 'both'. The hemisphere name. Used to construct the names of the surface file to be loaded. For 'both', see the information on the return value.
#'
#' @param force_hemilist logical, whether to return a \code{\link{hemilist}} even if the 'hemi' parameter is not set to 'both'
#'
#' @return the `fs.surface` instance, as returned by \code{\link[freesurferformats]{read.fs.surface}}. If parameter `hemi` is set to `both`, a named list with entries `lh` and `rh` is returned, and the values of are the respective surfaces. The mesh data structure used in `fs.surface` is a *face index set*.
#'
#' @family surface mesh functions
#'
#' @keywords internal
subject.surface <- function(subjects_dir, subject_id, surface = "white", hemi = "both", force_hemilist = FALSE) {

    if(!(hemi %in% c("lh", "rh", "both"))) {
        stop(sprintf("Parameter 'hemi' must be one of 'lh', 'rh' or 'both' but is '%s'.\n", hemi));
    }

    if(hemi == "both") {
        ret_list = list();
        ret_list$lh = subject.surface(subjects_dir, subject_id, surface, 'lh');
        ret_list$rh = subject.surface(subjects_dir, subject_id, surface, 'rh');
        return(ret_list);
    }

    surface_file = file.path(subjects_dir, subject_id, "surf", sprintf("%s.%s", hemi, surface));
    if(!file.exists(surface_file)) {
        stop(sprintf("Surface file '%s' for subject '%s' surface '%s' hemi '%s' cannot be accessed.\n", surface_file, subject_id, surface, hemi));
    }
    sf = freesurferformats::read.fs.surface(surface_file);

    if(force_hemilist) {
        return(hemilist.wrap(sf, hemi));
    } else {
        return(sf);
    }
}


#' @title Check whether object is a brainparc instance.
#'
#' @param x any `R` object
#'
#' @return logical, \code{TRUE} if parameter \code{x} is a brainparc instance (that is, has "brainparc" in its classes) and \code{FALSE} otherwise.
#'
#' @export
is.brainparc <- function(x) inherits(x, "brainparc")


#' @title Print description of a brainparce.
#'
#' @param x brainparc instance, with class `brainparc`.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @return Called for side effects.
#'
#' @export
print.brainparc <- function(x, ...) { # nocov start
    cat(sprintf("Brain parcellation for surface '%s'.\n * Left hemi with %d vertices and %d faces, right hemi with %d vertices and %d faces.\n", names(x$surfaces)[1], nrow(x$surfaces[[1]]$lh$vertices), nrow(x$surfaces[[1]]$lh$faces), nrow(x$surfaces[[1]]$rh$vertices), nrow(x$surfaces[[1]]$rh$faces)));
    num_atlases = length(names(x$annots));
    atl = "atlas";
    if(num_atlases > 1L) {
        atl = "atlases";
    }
    cat(sprintf(" * %d available brain %s: %s.\n", num_atlases, atl, paste(names(x$annots), collapse = ", ")));
    for (atlas in names(x$annots)) {
        cat(sprintf("  - atlas '%s' with %d lh regions and %d rh regions.\n", atlas, length(unique(x$annots[[atlas]]$lh)), length(unique(x$annots[[atlas]]$rh))));
    }
} # nocov end


#' @title Wrap data into a named hemi list.
#'
#' @param data something to wrap, typically some data for a hemisphere, e.g., a vector of morphometry data values. If NULL, the name will not be created.
#'
#' @param hemi character string, one of 'lh' or 'rh'. The name to use for the data in the returned list.
#'
#' @param hemilist optional \code{\link{hemilist}}, an existing hemilist to add the entry to. If left at the default value `NULL`, a new list will be created.
#'
#' @return a \code{\link{hemilist}}: a named list, with the 'data' in the name given by parameter 'hemi'
#'
#' @family hemilist functions
#'
#' @keywords internal
hemilist.wrap <- function(data, hemi, hemilist=NULL) {
    if(!(hemi %in% c("lh", "rh"))) {
        stop(sprintf("Parameter 'hemi' must be one of 'lh' or 'rh' but is '%s'.\n", hemi));
    }
    if(is.null(hemilist)) {
        ret_list = list();
    } else {
        ret_list = hemilist;
    }
    if(!is.null(data)) {
        ret_list[[hemi]] = data;
    }
    return(ret_list);
}

