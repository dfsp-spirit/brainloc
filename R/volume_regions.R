

#' @title Compute the center of mass of the given points.
#'
#' @param coords numeric \code{nx3} matrix if point x,y,z coordinates. A single coordinate can be passed as a vector of length 3, and will be returned as is.
#'
#' @param weights numerical vector of length \code{n}, weights for the points in \code{coords}. Assumed to be all \code{1.0} if omitted.
#'
#' @return vector of length 3, the center of mass.
#'
#' @keywords internal
center_of_mass <- function(coords, weights=NULL) {
    if(is.vector(coords)) {
        if(length(coords) == 3L) {
            return(coords);
        } else {
            stop("If 'coords' is a vector, it must have length exactly 3.");
        }
    }
    if(ncol(coords) != 3L) {
        stop("Parameter 'coords' must have exactly 3 columns");
    }
    nc = nrow(coords); # nc = number of coords
    if(is.null(weights)) {
        weights = rep(1.0, nc);
    } else {
        if(length(weights) != nc) {
            stop(sprintf("Received %d coordinates but %d weights. Counts must match.\n", nc, length(weights)));
        }
    }
    mass = sum(weights);
    x = sum(coords[,1] * weights);
    y = sum(coords[,2] * weights);
    z = sum(coords[,3] * weights);
    return(c(x/mass, y/mass, z/mass));
}


#' @title Get names of regions in labeled volume from color lookup table (colorLUT).
#'
#' @description Given a brain volume with integer labels as voxel values and a suitable colorLUT, extract the region names as a named list.
#'
#' @param volatlas 3d integer array, the brain volume segmentation. Can be an \code{fs.volume} instance. Typically something like \code{'subject/mri/aseg.mgz'}.
#'
#' @param colorlut a color lookup table data.frame, as returned by \code{\link[freesurferformats]{read.fs.colortable}}. Can be \code{NULL} if you do not have any, in that case the region names will just be 'region_x', where 'x' is the integer region code in the volatlas.
#'
#' @param ignore_not_in_lut logical, whether to ignore regions that are not found in the colorLUT. If TRUE, the regions will not appear in the named output list. If FALSE, they will appear as 'region_x', where 'x' is the integer region code in the volatlas.
#'
#' @param warn_not_in_lut logical, whether to print a warning if some regions of the volume do not occur in the colorLUT.
#'
#' @return named list, the keys are character strings: the region names from the colorLUT, and the values are integers: the region codes from the volume.
#'
#' @export
name_regions <- function(volatlas, colorlut, ignore_not_in_lut=FALSE, warn_not_in_lut=FALSE) {
    if(is.character(volatlas)) {
        volatlas = freesurferformats::read.fs.volume(volatlas);
    }
    if(freesurferformats::is.fs.volume(volatlas)) {
        volatlas = volatlas$data;
    }
    if(is.character(colorlut)) {
        colorlut = freesurferformats::read.fs.colortable(colorlut);
    }
    if(! is.array(volatlas)) {
        stop("Parameter volatlas must be a 3D array of integers.");
    }
    num_dim = length(dim(drop(volatlas)));
    if(num_dim != 3L) {
        stop("Volume in parameter 'volatlas' must be 3D (or 4D with a single frame).");
    }
    region_codes = unique(as.integer(volatlas));
    named_reg = list();
    codes_not_found = c();
    for(reg_code in region_codes) {
        found = FALSE;
        if(! is.null(colorlut)) {
            if(reg_code %in% colorlut$struct_index) {
                reg_name = colorlut[colorlut$struct_index == reg_code, ]$struct_name;
                named_reg[reg_name] = reg_code;
                found = TRUE;
            }
        }

        if(! found) {
            codes_not_found = c(codes_not_found, reg_code);
            if(! ignore_not_in_lut) {
                reg_name = sprintf("region_%d", reg_code);
                named_reg[reg_name] = reg_code;
            }
        }
    }
    if(length(codes_not_found) > 0L) {
        if(warn_not_in_lut) {
            warning(sprintf("There were %d region codes in the segmentation that are not listed in the colorLUT: %s.\n", length(codes_not_found), paste(codes_not_found, collapse = ", ")));
        }
    }
    return(named_reg);
}


#' @title Find center of mass of regions in a volume segmentation.
#'
#' @inheritParams name_regions
#'
#' @param named_regions the regions to consider, a named list where the keys are character strings: the region names from the colorLUT, and the values are integers: the region codes from the volume. See \code{\link{name_regions}} to get one. Can be \code{NULL}, in which case all regions in the volatlas will be used, and the region names will appear as 'region_x', where 'x' is the integer region code in the volatlas.
#'
#' @param regions
#'
#' @examples
#' \dontrun{
#' fsh = brainloc:::fs.home();
#' lutfile = file.path(fsh, "FreeSurferColorLUT.txt");
#' segfile = file.path(fsh, 'subjects', 'fsaverage', 'mri', 'aseg.mgz');
#' named_regions = name_regions(segfile, lutfile);
#' segmentation_centers(segfile, named_regions);
#' }
#'
#' @export
segmentation_centers <- function(volatlas, named_regions=NULL, vox2ras=diag(4)) {
    if(is.character(volatlas)) {
        volatlas = freesurferformats::read.fs.volume(volatlas);
    }
    if(freesurferformats::is.fs.volume(volatlas)) {
        if(is.null(vox2ras)) {
            vox2ras = freesurferformats::mghheader.vox2ras(volatlas$header);
        }
        volatlas = volatlas$data;
    }
    if(is.null(named_regions)) {
        named_regions = named_regions(volatlas, colorlut = NULL);
    }
    if(is.null(vox2ras)) {
        vox2ras=diag(4);
    }

    for(reg_name in names(named_regions)) {
        reg_code = named_regions[[reg_name]];
        cat(sprintf("Handling region '%s' with code %d.\n", reg_name, reg_code));
    }
}

