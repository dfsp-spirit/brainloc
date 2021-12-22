

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
