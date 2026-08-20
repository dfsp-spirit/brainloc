# Find closest regions to coordinate using Euclidean or geodesic distance.

Finds the closest atlas regions according to the brain surface
parcellation 'brainparc' for the given query vertex coordinates. The
coordinates must be in the surface space of the surface included in the
'brainparc'.

## Usage

``` r
coord_closest_regions(
  brainparc,
  coordinate,
  linkage = "single",
  distance = "euclidean",
  silent = getOption("brainloc.silent", default = FALSE),
  num_regions_to_report = 5L
)
```

## Arguments

- brainparc:

  a brain parcellation, see functions like
  [`brainparc_fs`](https://dfsp-spirit.github.io/brainloc/reference/brainparc_fs.md)
  to get one.

- coordinate:

  `nx3` numerical matrix or vector of length 3, the query point
  coordinates.

- linkage:

  character string, one of `"single"` or `"centroid"`. Defines how the
  distance from a vertex to a region of vertices is computed.
  `"single"`: Euclidean distance from query vertex to the closest vertex
  of the atlas region. `"centroid"`: Euclidean distance from query
  vertex to the mean of the vertex coordinates of the atlas region.

- distance:

  character string, one of `"euclidean"` or `"geodesic"`. The latter is
  only supported with `linkage = 'single'`.

- silent:

  logical, whether to suppress console messages.

- num_regions_to_report:

  scalar integer, the number of closest regions to report (will be
  reported in increasing order, by distance). Pass `NULL` to get all
  regions.

## Note

This is a wrapper around
[`vertex_closest_regions`](https://dfsp-spirit.github.io/brainloc/reference/vertex_closest_regions.md).
It first computes the surface vertex closest to the given coordinate
(using
[`coord_closest_vertex`](https://dfsp-spirit.github.io/brainloc/reference/coord_closest_vertex.md)),
then runs
[`vertex_closest_regions`](https://dfsp-spirit.github.io/brainloc/reference/vertex_closest_regions.md)
with that vertices coordinates.

## See also

[`vertex_closest_regions`](https://dfsp-spirit.github.io/brainloc/reference/vertex_closest_regions.md)
is faster if you have a vertex index for the surface in 'brainparc'
instead of a coordinate.

## Examples

``` r
if (FALSE) { # \dontrun{
bp = brainparc_fs(fsbrain::fsaverage.path(), "fsaverage", atlas="aparc");
coord_closest_regions(bp, bp$surfaces$white$lh$vertices[1:3,]);
} # }

```
