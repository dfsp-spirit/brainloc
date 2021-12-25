


#' @title Convert fs.surface to tmesh3d using freesurferformats of fsbrain.
#'
#' @param surface an fs.surface instance.
#'
#' @return \code{tmesh} instance, as used in the \code{rgl} package.
#'
#' @keywords internal
fs.surface.to.tmesh3d <- function(surface) {
    if( ! freesurferformats::is.fs.surface(surface)) {
        stop("Parameter 'surface' must be an instance of freesurferformats::fs.surface.");
    }
    tmesh = list("material"=list(), "normals"=NULL, "texcoords"=NULL, "meshColor"="vertices");
    class(tmesh) = c("mesh3d", "shape3d");
    tmesh$vb = t(cbind(surface$vertices, 1L)); # Transform vertex coords to homogeneous and swap rows/columns
    tmesh$it = t(surface$faces); # swap only
    return(tmesh);
}


#' @title Return FreeSurfer path.
#'
#' @return the FreeSurfer path, typically what the environment variable `FREESURFER_HOME` points to.
#'
#' @note This function will stop (i.e., raise an error) if the directory cannot be found. It calls \code{\link{find.freesurferhome}} internally, see there for more details.
#'
#' @family FreeSurfer helper functions
#'
#' @export
fs_home <- function() {
    return(find.freesurferhome(mustWork=TRUE));
}


#' @title Find subjects_dir or stop.
#'
#' @description Try a number of ways to find the subjects dir, using environment variables and knowledge on common installation paths of FreeSurfer.
#'
#' @param mustWork whether to stop if no subjects_dir can be found.
#'
#' @param allow_download logical, whether to allow downloading the data in case it is not found locally.
#'
#' @inheritParams sjd_demo
#'
#' @return character string, the subjects directory. If mustWork is FALSE, it will return NULL if no directory was found. If mustWork is TRUE and nothing is found, it stops.
#'
#' @family FreeSurfer helper functions
#'
#' @export
get_subjects_dir <- function(mustWork = TRUE, allow_download=FALSE, accept_freesurfer_license = FALSE) {

    guessed_path = getOption("brainloc.subjects_dir");
    if(! is.null(guessed_path)) {
        if(dir.exists(guessed_path)) {
            return(guessed_path);
        } else {
            stop(sprintf("The path given by 'getOption('brainloc.subjects_dir')', '%s', does not exist. Please fix or unset.\n", guessed_path));
        }
    }

    env_subjects_dir=Sys.getenv("SUBJECTS_DIR");
    if(nchar(env_subjects_dir) > 0) {
        guessed_path = file.path(env_subjects_dir);
        if(dir.exists(guessed_path)) {
            return(guessed_path);
        }
    }

    # If not found, try based on fs home.
    fs_home = find.freesurferhome(mustWork = FALSE);
    if(fs_home$found) {
        fs_home = fs_home$found_at;
        guessed_path = file.path(fs_home, 'subjects');
        if(dir.exists(guessed_path)) {
            return(guessed_path);
        }
    }

    if(allow_download) {
        if(! accept_freesurfer_license) {
            stop("You must read and accept the FreeSurfer license, and then set parameter 'accept_freesurfer_license' to TRUE to download the FreeSurfer template data.");
        }
        return(sjd_demo(accept_freesurfer_license = TRUE));
    }

    if(mustWork) {
        stop("Could not find subjects_dir, please set environment variable SUBJECTS_DIR.");
    } else {
        return(NULL);
    }
}

#' @title Check whether a FreeSurfer installation can be found on the system.
#'
#' @return logical, whether a FreeSurfer installation was found.
#'
#' @family FreeSurfer helper functions
#'
#' @export
has_fs <- function() {
    return(find.freesurferhome(mustWork = FALSE)$found);
}

#' @title Find the FREESURFER_HOME directory on disk.
#'
#' @description Try to find directory containing the FreeSurfer installation, based on environment variables and *educated guessing*.
#'
#' @param mustWork logical. Whether the function should with an error stop if the directory cannot be found. If this is TRUE, the return value will be only the 'found_at' entry of the list (i.e., only the path of the FreeSurfer installation dir).
#'
#' @return named list with the following entries: "found": logical, whether it was found. "found_at": Only set if found=TRUE, the path to the FreeSurfer installation directory (including the directory itself). See 'mustWork' for important information.
#'
#' @seealso \code{\link{fs_home}}, \code{\link{has_fs}}
#'
#' @family FreeSurfer helper functions
#'
#' @keywords internal
find.freesurferhome <- function(mustWork=FALSE) {
    ret = list();
    ret$found = FALSE;

    fs_home=Sys.getenv("FREESURFER_HOME");
    if(nchar(fs_home) > 0) {
        guessed_path = file.path(fs_home);
        if(dir.exists(guessed_path)) {
            ret$found = TRUE;
            ret$found_at = guessed_path;
        }
    }

    # Check in some typical paths
    if(! ret$found) {
        if(tolower(Sys.info()[["sysname"]]) == 'darwin') {
            search_paths = c("/Applications/freesurfer");
        } else if(tolower(Sys.info()[["sysname"]]) == 'linux') {
            search_paths = c("/usr/local/freesurfer", "/opt/freesurfer");
        } else {
            # Windows, needed for AppVeyor
            search_paths = c();
        }

        user_home = Sys.getenv("HOME");
        if(nchar(user_home) > 0) {
            search_paths = c(search_paths, file.path(user_home, 'freesurfer'), file.path(user_home, 'software', 'freesurfer'), file.path(user_home, 'opt', 'freesurfer'));
        }

        for(sp in search_paths) {
            if(dir.exists(sp)) {
                ret$found = TRUE;
                ret$found_at = sp;
            }
        }

    }

    if(mustWork) {
        if(ret$found) {
            return(ret$found_at);
        } else {
            stop(sprintf("Could not find FreeSurfer installation dir and parameter 'mustWork' is TRUE. Please set the environment variables by installing and configuring FreeSurfer.\n"));
        }
    }

    return(ret);
}


#' @title Create hemilist.
#'
#' @description A hemilist is a named list which has at least one of the following fields: \code{lh} and \code{rh}. The fields can store arbitrary values.
#'
#' @param lh_data any R object, the lh value.
#'
#' @param rh_data any R object, the rh value.
#'
#' @return a hemilist
#'
#' @note This function is trivial and mainly exists to document what a hemilist is (so this explanation can be linked to from the documentation of other functions which use the concept).
#'
#' @family hemilist functions
#'
#' @export
hemilist <- function(lh_data=NULL, rh_data=NULL) {
    return(list('lh' = lh_data, 'rh' = rh_data));
}


#' @title Find out whether x is a hemilist
#'
#' @param x any R object
#'
#' @return logical, whether x is a hemilist
#'
#' @family hemilist functions
#'
#' @keywords internal
is.hemilist <- function(x) {
    if(! is.list(x)) {
        return(FALSE);
    }
    if(!("lh" %in% names(x) | "rh" %in% names(x))) {
        return(FALSE)
    }
    return(TRUE);
}



