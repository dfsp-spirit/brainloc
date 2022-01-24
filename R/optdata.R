#' @title Download the FreeSurfer v6 fsaverage subject.
#'
#' @description Download some relevant files from the FreeSurfer v6 fsaverage subject. The files are subject to the FreeSurfer software license, see parameter 'accept_freesurfer_license' for details. This data is not required for the package to work. If you are working on a machine that has FreeSurfer installed, you already have this data anyways and do not need to download it. If not, it is very convenient to have it if you are using the fsaverage template subject to analyze your standard space data, as it is required for visualization of such data.
#'
#' @param accept_freesurfer_license logical, whether you accept the FreeSurfer license for fsaverage, available at https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense. Defaults to FALSE.
#'
#' @return Named list. The list has entries: "available": vector of strings. The names of the files that are available in the local file cache. You can access them using get_optional_data_file(). "missing": vector of strings. The names of the files that this function was unable to retrieve.
#'
#' @keywords internal
download_fsaverage <- function(accept_freesurfer_license=FALSE) {

    if(! accept_freesurfer_license) {
        cat(sprintf("Nothing downloaded. You have to accept the FreeSurfer license to download and use fsaverage.\n"));
        cat(sprintf("Read the license at https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense and set parameter 'accept_freesurfer_license' to TRUE if you accept it.\n"));
        return(invisible(NULL));
    }

    pkg_info = pkgfilecache::get_pkg_info("brainloc");
    base_path_fsaverage = c('subjects_dir', 'fsaverage');
    local_filenames = list(c(base_path_fsaverage, 'label', 'lh.aparc.a2009s.annot'),
                           c(base_path_fsaverage, 'label', 'rh.aparc.a2009s.annot'),
                           c(base_path_fsaverage, 'label', 'lh.aparc.annot'),
                           c(base_path_fsaverage, 'label', 'rh.aparc.annot'),
                           c(base_path_fsaverage, 'label', 'lh.cortex.label'),
                           c(base_path_fsaverage, 'label', 'rh.cortex.label'),
                           c(base_path_fsaverage, 'mri', 'brain.mgz'),
                           c(base_path_fsaverage, 'surf', 'lh.white'),
                           c(base_path_fsaverage, 'surf', 'rh.white'),
                           c(base_path_fsaverage, 'surf', 'lh.pial'),
                           c(base_path_fsaverage, 'surf', 'rh.pial'),
                           c(base_path_fsaverage, 'surf', 'lh.inflated'),
                           c(base_path_fsaverage, 'surf', 'rh.inflated'),
                           c(base_path_fsaverage, 'surf', 'lh.curv'),
                           c(base_path_fsaverage, 'surf', 'rh.curv'),
                           c(base_path_fsaverage, 'ext', 'FreeSurferColorLUT.txt'),
                           c(base_path_fsaverage, 'LICENSE')
    );



    md5sums = c('b4310b1e4435defaf27fc7ee98199e6a',
                '6077dc6cb42dd8c48bb382672d65743c',
                'bf0b488994657435cdddac5f107d21e8',
                '8f504caddedfde367a40501da6222809',
                '578f81e9946a76eb1c42d897d07da4a7',
                'c8f59de23e9f90f18e96e9d037e42799',
                'b8bc4b5854f2d5e66d5c4f95d4f9cf63',
                'cbffce8198e0e10c17f79f6ae0454af5', # lh.white
                '1159a9ee160b1b0c76e0bb9ae789b9be',
                'c53c1f70ae8971e1c04bd19e3277fa14',
                '71f11c33db672360d7589c7dbd0e4a3f',
                '95df985980d7eefa009ac104589ee3c5',
                'bb4d58289aefcdf8d017e45e531c4807',
                '3e81598a5ac0546443ec37d0ac477c80',
                '76ad91d2488de081392313ad5a87fafb',
                'a3735566ef949bd4d7ed303837cc5e77',  # color LUT
                'b39610adfe02fdce2ad9d30797c567b3'    # LICENSE
    );



    ext_url_subject_part_fsaverage = 'subjects_dir/fsaverage/';
    ext_url_parts_each_subject = c('label/lh.aparc.a2009s.annot',
                                   'label/rh.aparc.a2009s.annot',
                                   'label/lh.aparc.annot',
                                   'label/rh.aparc.annot',
                                   'label/lh.cortex.label',
                                   'label/rh.cortex.label',
                                   'mri/brain.mgz',
                                   'surf/lh.white',
                                   'surf/rh.white',
                                   'surf/lh.pial',
                                   'surf/rh.pial',
                                   'surf/lh.inflated',
                                   'surf/rh.inflated',
                                   'surf/lh.curv',
                                   'surf/rh.curv',
                                   'ext/FreeSurferColorLUT.txt',
                                   'LICENSE'
    );
    ext_urls = paste(ext_url_subject_part_fsaverage, ext_url_parts_each_subject, sep='');
    base_url = 'http://rcmd.org/projects/nitestdata/';
    urls = paste(base_url, ext_urls, sep='');

    cfiles = pkgfilecache::ensure_files_available(pkg_info, local_filenames, urls, md5sums=md5sums);
    cfiles$file_status = NULL; # not exposed to end user
    return(invisible(cfiles));
}


#' @title Download the FreeSurfer v6 low-resolution fsaverage3 subject.
#'
#' @description Download some relevant files from the FreeSurfer v6 fsaverage3 subject. The files are subject to the FreeSurfer software license, see parameter 'accept_freesurfer_license' for details. This data is not required for the package to work. If you are working on a machine that has FreeSurfer installed, you already have this data anyways and do not need to download it. Also downloads data for subject1 that has been mapped to fsaverage.
#'
#' @inheritParams download_fsaverage
#'
#' @return Named list. The list has entries: "available": vector of strings. The names of the files that are available in the local file cache. You can access them using get_optional_data_file(). "missing": vector of strings. The names of the files that this function was unable to retrieve.
#'
#' @note The subject fsaverage3 is a downsampled (low mesh resolution) version of the standard fsaverage. If you never heard about fsaverage3, you do not need it. You will have to manually re-sample your data in FreeSurfer if you want to use it with fsaverage3.
#'
#' @keywords internal
download_fsaverage3 <- function(accept_freesurfer_license=FALSE) {

    if(! accept_freesurfer_license) {
        cat(sprintf("Nothing downloaded. You have to accept the FreeSurfer license to download and use fsaverage.\n"));
        cat(sprintf("Read the license at https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense and set parameter 'accept_freesurfer_license' to TRUE if you accept it.\n"));
        return(invisible(NULL));
    }

    pkg_info = pkgfilecache::get_pkg_info("brainloc");
    base_path_fsaverage3 = c('subjects_dir', 'fsaverage3');
    base_path_subject1 = c('subjects_dir', 'subject1');
    local_filenames = list(c(base_path_fsaverage3, 'label', 'lh.cortex.label'),
                           c(base_path_fsaverage3, 'label', 'rh.cortex.label'),
                           c(base_path_fsaverage3, 'surf', 'lh.white'),
                           c(base_path_fsaverage3, 'surf', 'rh.white'),
                           c(base_path_fsaverage3, 'LICENSE'),
                           c(base_path_subject1, 'surf', 'lh.thickness.fwhm0.fsaverage3.mgz'),
                           c(base_path_subject1, 'surf', 'rh.thickness.fwhm0.fsaverage3.mgz')
    );



    md5sums = c('49a367e65ec7ecffbb721404b274fb3f', # fsaverage3 lh.cortex
                '76e2d42894351427405cc01ab351719b',
                'b014033974bc5b4deb8b54dc140abda8',
                '09a133fd8499f3192e051bdbd8bec6e8',
                'b39610adfe02fdce2ad9d30797c567b3',    # LICENSE fsaverage3
                'd191f6833d1d36016b30504fed1ce138',  # subject1 lh.thickness.fwhm0.fsaverage3.mgz
                'e874f8dc149fd11842f117c300d1a964'   # subject1 rh.thickness.fwhm0.fsaverage3.mgz
    );



    ext_url_subject_part_fsaverage3 = 'subjects_dir/fsaverage3/';
    ext_url_parts_fsaverage3 = c('label/lh.cortex.label',
                                 'label/rh.cortex.label',
                                 'surf/lh.white',
                                 'surf/rh.white',
                                 'LICENSE'
    );
    ext_urls_fsaverage3 = paste(ext_url_subject_part_fsaverage3, ext_url_parts_fsaverage3, sep='');

    ext_url_subject_part_subject1 = 'subjects_dir/subject1/';
    ext_url_parts_subject1 = c('surf/lh.thickness.fwhm0.fsaverage3.mgz',
                               'surf/rh.thickness.fwhm0.fsaverage3.mgz'
    );
    ext_urls_subject1 = paste(ext_url_subject_part_subject1, ext_url_parts_subject1, sep='');
    ext_urls = c(ext_urls_fsaverage3, ext_urls_subject1);
    base_url = 'http://rcmd.org/projects/nitestdata/';
    urls = paste(base_url, ext_urls, sep='');

    cfiles = pkgfilecache::ensure_files_available(pkg_info, local_filenames, urls, md5sums=md5sums);
    cfiles$file_status = NULL; # not exposed to end user
    return(invisible(cfiles));
}


#' @title Download the talairach volume and labels from talairach.org.
#'
#' @description Download the talairach volume \code{talairach.nii} and labels \code{labels.txt} from talairach.org.
#'
#' @return Named list. The list has entries: "available": vector of strings. The names of the files that are available in the local file cache. You can access them using get_optional_data_file(). "missing": vector of strings. The names of the files that this function was unable to retrieve.
#'
#' @note This function requires and internet connection. Files will only be downloaded if they are not already available on the local computer.
#'
#' @keywords internal
download_talairach <- function(accept_talairach_usage=FALSE) {

    if(! accept_talairach_usage) {
        cat(sprintf("Nothing downloaded. You have to accept the talairach.org usage conditions to download the data.\n"));
        cat(sprintf("Read which publications to cite when using these data at http://talairach.org/ and set parameter 'accept_talairach_usage' to TRUE if you accept it.\n"));
        return(invisible(NULL));
    }

    pkg_info = pkgfilecache::get_pkg_info("brainloc");
    base_path_tal = c('talairach');
    local_filenames = list(c(base_path_tal, 'talairach.nii'),
                           c(base_path_tal, 'labels.txt'));

    md5sums = c('6ddca3a02fd8a90e1b6daf827d4a9d98', # talairach.nii
                'f4c0bc758696ea52be6c11848df4786d'  # labels.txt
    );

    ext_urls = c('talairach.nii', 'labels.txt'); # the files are directly in the main dir on the web server, no subdirs/other URL parts needed
    base_url = 'http://www.talairach.org/';
    urls = paste(base_url, ext_urls, sep='');

    cfiles = pkgfilecache::ensure_files_available(pkg_info, local_filenames, urls, md5sums=md5sums);
    cfiles$file_status = NULL; # not exposed to end user
    return(invisible(cfiles));
}



#' @title Get file names available in package cache.
#'
#' @description Get file names of optional data files which are available in the local package cache. You can access these files with get_optional_data_file().
#'
#' @return vector of strings. The file names available, relative to the package cache.
#'
#' @keywords internal
list_optional_data <- function() {
    pkg_info = pkgfilecache::get_pkg_info("brainloc");
    return(pkgfilecache::list_available(pkg_info));
}


#' @title Access a single file from the package cache by its file name.
#'
#' @param filename, string. The filename of the file in the package cache.
#'
#' @param mustWork, logical. Whether an error should be created if the file does not exist. If mustWork=FALSE and the file does not exist, the empty string is returned.
#'
#' @return string. The full path to the file in the package cache or the empty string if there is no such file available. Use this in your application code to open the file.
#'
#' @keywords internal
get_optional_data_filepath <- function(filename, mustWork=TRUE) {
    pkg_info = pkgfilecache::get_pkg_info("brainloc");
    return(pkgfilecache::get_filepath(pkg_info, filename, mustWork=mustWork));
}


#' @title Delete all data in the package cache.
#'
#' @return integer. The return value of the unlink() call: 0 for success, 1 for failure. See the unlink() documentation for details.
#'
#' @keywords internal
delete_all_optional_data <- function() {
    pkg_info = pkgfilecache::get_pkg_info("brainloc");
    return(pkgfilecache::erase_file_cache(pkg_info));
}


#' @title Download FreeSurfer template and demo data if needed and return its path (subjects_dir).
#'
#' @description This is a wrapper around \code{download_fsaverage()} and \code{download_fsaverage3())}. It will download the data from the internet unless it already exists locally.
#'
#' @param accept_freesurfer_license logical, whether you want to also download fsaverage and fsaverage3, and accept the FreeSurfer license for fsaverage and fsaverage3, available at https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense. Defaults to FALSE. If FALSE, nothing will be downloaded.
#'
#' @return character string, the path to the 'subjects_dir' directory within the downloaded template data directory.
#'
#' @note This function will stop if the data cannot be accessed, i.e., the 'subjects_dir' does not exist after trying to download the data.
#'
#' @keywords internal
sjd_demo <- function(accept_freesurfer_license=FALSE) {
    download_fsaverage(accept_freesurfer_license);
    download_fsaverage3(accept_freesurfer_license);
    return(get_optional_data_filepath("subjects_dir"));
}




