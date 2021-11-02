

#' @title Create a brain parcellation from a FreeSurfer annotation.
#'
#' @inheritParams subject.surface
#'
#' @inheritParams subject.annot
#'
#' @export
brainparc_fs <- function(subjects_dir, subject_id, surface="white", atlas=c("aparc", "aparc.a2009s")) {
    surfaces = subject.surface(subjects_dir, subject_id, hemi="both", surface = surface);
    ret = list("surfaces"=list(), "annots"=list());
    ret$surfaces[[surface]] = surfaces;
    for(atlas_name in atlas) {
        ret$annots[[atlas_name]] = subject_annot_simple_fs(subjects_dir, subject_id, atlas=atlas_name);
    }
    class(ret) = c(class(ret), "brainparc");
    return(ret);
}

#' @keywords internal
get_surface <- function(brainparc) {
    return(brainparc$surfaces[[1]]);
}


#' @title Return simplified annotation: the region names only.
#'
#' @inheritParams subject.annot
#'
#' @keywords internal
subject_annot_simple_fs <- function(subjects_dir, subject_id, atlas) {
    lh_annot = subject.annot(subjects_dir, subject_id, hemi="lh", atlas);
    rh_annot = subject.annot(subjects_dir, subject_id, hemi="rh", atlas);
    return(list("lh"=lh_annot$label_names, "rh"=rh_annot$label_names));
}



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
#' @param force_hemilist logical, whether to return a hemilist even if the 'hemi' parameter is not set to 'both'
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


#' @title Wrap data into a named hemi list.
#'
#' @param data something to wrap, typically some data for a hemisphere, e.g., a vector of morphometry data values. If NULL, the name will not be created.
#'
#' @param hemi character string, one of 'lh' or 'rh'. The name to use for the data in the returned list.
#'
#' @param hemilist optional hemilist, an existing hemilist to add the entry to. If left at the default value `NULL`, a new list will be created.
#'
#' @return a hemilist: a named list, with the 'data' in the name given by parameter 'hemi'
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

