# Find closest surface vertex to a point, using Euclidean distance.

Find closest surface vertex to a point, using Euclidean distance.

## Usage

``` r
coord_closest_vertex(coordinate, surfaces)
```

## Arguments

- coordinate:

  `nx3` numerical matrix or vector of length 3, the query point
  coordinates.

- surfaces:

  [`hemilist`](https://dfsp-spirit.github.io/brainloc/reference/hemilist.md)
  of `fs.surface` instances

## Value

a data.frame with columns named 'query_x', 'query_y', 'query_z',
'lh_closest_vertex', 'lh_distance', 'rh_closest_vertex', 'rh_distance',
'both_closest_vertex', 'both_distance', 'both_hemi'.

## Examples

``` r
if (FALSE) { # \dontrun{
bp = brainparc_fs(fsbrain::fsaverage.path(), "fsaverage", atlas="aparc");
query_coords = matrix(seq.int(9), ncol = 3, byrow = TRUE);
coord_closest_vertex(query_coords, bp$surfaces$white);
} # }
```
