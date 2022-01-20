# Functions for querying Talairach labels for coordinates in Talairach space.
# These functions require the Talairach lookup data (volume with label indices, and mapping of indices to label strings) from talairach.org.
# All

#' @title Retrieve Talairach labels for Talairach coordinates from \code{talairach.org} files.
#'
#' @param talfile the NIFTI Talairach label volume file from \code{talairach.org}. It is named \code{talairach.nii} on the website. You can also gzip it and use a much smaller \code{.nii.gz} version of the file.
#'
#' @param lookup_table_file the Talairach label index text file from \code{talairach.org}. It is named \code{labels.txt} on the website.
#'
#' @param tal_coords \code{nx3} numeric matrix, the \code{n} query Talairach coordinates for which you want to retrieve the labels.
#'
#' @return data.frame describing labels for the coordinates. The following columns are included: \code{cx,cy,cz}: the query coordinates, as given in parameter 'tal_coords'. \code{v1,v2,v3}: The voxel indices in the talairach.nii file that map to the query coordinates. \code{label_lvl1,...,label_lvl5}: the label strings for the location, in 5 levels (parsed from \code{label_full} for you as a convenience). \code{label_full}: The raw, full label string for the location.
#'
#' @note You need to download the required files from \code{talairach.org}. It can be found directly on the \code{home} page in the section \code{Talairach Label Data}.
#'
#' @note When using this function or the Talairach.org data in your research, please cite the following two publications: \code{Lancaster JL, Woldorff MG, Parsons LM, Liotti M, Freitas CS, Rainey L, Kochunov PV, Nickerson D, Mikiten SA, Fox PT, "Automated Talairach Atlas labels for functional brain mapping". Human Brain Mapping 10:120-131, 2000.} and \code{Lancaster JL, Rainey LH, Summerlin JL, Freitas CS, Fox PT, Evans AC, Toga AW, Mazziotta JC. Automated labeling of the human brain: A preliminary report on the development and evaluation of a forward-transform method. Hum Brain Mapp 5, 238-242, 1997.}
#'
#' @keywords internal
get_talairach_label <- function(talfile="~/develop/talairach/talairach.nii", lookup_table_file="~/develop/talairach/labels.txt", tal_coords=matrix(seq.int(15), nrow=5, ncol=3)) {

    if(! file.exists(talfile)) {
        stop(sprintf("Please download the 'talairach.nii' file from talairach.org and specify the correct path. Expected file not found at '%s'.\n", talfile));
    }
    if(! file.exists(lookup_table_file)) {
        stop(sprintf("Please download the 'labels.txt' file from talairach.org and specify the correct path. Expected file not found at '%s'.\n", lookup_table_file));
    }

    if(is.vector(tal_coords)) {
        tal_coords = matrix(tal_coords, ncol = 3L);
    }


    tal = freesurferformats::read.fs.volume(talfile, with_header = TRUE);
    lab = read.table(lookup_table_file, sep = "\t", col.names = c("index", "region"));

    taldata = drop(tal$data);

    r2v = freesurferformats::mghheader.ras2vox(tal);
    voxels = floor(freesurferformats::doapply.transform.mtx(tal_coords, r2v, as_mat = TRUE));

    label_indices = taldata[voxels];
    voxel_labels = lab$region[label_indices]; # each label is single string that looks like: ''
    voxel_labels_split = strsplit(voxel_labels, split = ".", fixed = TRUE);

    res = data.frame("cx"=tal_coords[,1], "cy"=tal_coords[,2], "cz"=tal_coords[,3], "v1"=voxels[,1], "v2"=voxels[,2], "v3"=voxels[,3], stringsAsFactors = FALSE);

    ncoords = length(voxel_labels_split); # Could also use the number of query voxels.
    nlevels = 5L;
    for(level_idx in seq(nlevels)) {
        level_names = c();
        for(coord_idx in seq(ncoords)) {
            level_names = c(level_names, voxel_labels_split[[coord_idx]][level_idx]);
            key = paste("label_lvl", level_idx, sep = "");
        }
        res[[key]] = level_names;
    }
    res$label_full = voxel_labels;
    return(res);
}

