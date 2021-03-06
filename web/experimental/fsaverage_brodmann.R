#!/usr/bin/env Rscript


library("freesurferformats");
library("fsbrain");
library("brainloc");

#' @title Try to create a Broadmann area parcellation for the FreeSurfer fsaverage subject.
#'
#' @description This bold demo function tries to create a visualization of Brodmann areas on the FreeSurfer MNI305 space fsaverage subject. As the Brodmann data comes from a NIFTI volume in Talairach space and goes via MNI152 (piecewise-linear) to MNI305, the result is expected to be very rough. A better approach would most likely be to map from Talairach to MIN152 volume, then MNI305 volume, the project to the surface.
#'
#' @seealso \code{https://surfer.nmr.mgh.harvard.edu/fswiki/BrodmannAreaMaps}
#'
#' @return hemilist of level 5 Talairach label for the fsaverage vertices, but it is rather called for the side effect of plotting the atlas.
fsaverage_brodmann_atlas <- function() {
    sjd = fsbrain::fsaverage.path(T);
    sj = "fsaverage";

    surfaces = fsbrain::subject.surface(sjd, sj, surface="pial"); # Note that which surface is used can make a large difference.
    talairach_coords = list("lh"= brainloc::coord_MNI305_info(surfaces$lh$vertices)$talairach, "rh"= brainloc::coord_MNI305_info(surfaces$rh$vertices)$talairach);
    talairach_labels = list("lh"=brainloc::get_talairach_label(talairach_coords$lh), "rh"=brainloc::get_talairach_label(talairach_coords$rh));
    labels_int = list("lh"=as.integer(as.factor(talairach_labels$lh$label_lvl5)), "rh"=as.integer(as.factor(talairach_labels$rh$label_lvl5)) );
    fsbrain::vis.data.on.fsaverage(morph_data_lh = labels_int$lh, morph_data_rh = labels_int$rh, views = "si"); # show the atlas.
    return(talairach_labels);
}

tl = fsaverage_brodmann_atlas();

