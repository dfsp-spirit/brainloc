# Compute the center of mass of the given points.

Compute the center of mass of the given points.

## Usage

``` r
center_of_mass(coords, weights = NULL)
```

## Arguments

- coords:

  numeric `nx3` matrix if point x,y,z coordinates. A single coordinate
  can be passed as a vector of length 3, and will be returned as is.

- weights:

  numerical vector of length `n`, weights for the points in `coords`.
  Assumed to be all `1.0` if omitted.

## Value

vector of length 3, the center of mass.
