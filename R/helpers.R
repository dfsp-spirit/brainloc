


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
#' @keywords internal
fs.home <- function() {
    return(find.freesurferhome(mustWork=TRUE));
}


#' @title Find the FREESURFER_HOME directory on disk.
#'
#' @description Try to find directory containing the FreeSurfer installation, based on environment variables and *educated guessing*.
#'
#' @param mustWork logical. Whether the function should with an error stop if the directory cannot be found. If this is TRUE, the return value will be only the 'found_at' entry of the list (i.e., only the path of the FreeSurfer installation dir).
#'
#' @return named list with the following entries: "found": logical, whether it was found. "found_at": Only set if found=TRUE, the path to the FreeSurfer installation directory (including the directory itself). See 'mustWork' for important information.
#'
#' @seealso \code{\link{fs.home}}
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

