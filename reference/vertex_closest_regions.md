# Find closest regions to vertex using Euclidean or geodesic distance.

Finds the closest atlas regions according to the brain surface
parcellation 'brainparc' for the given query vertices.

## Usage

``` r
vertex_closest_regions(
  brainparc,
  vertices,
  hemis,
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

- vertices:

  integer vector, the query vertex indices (from the surface in the
  `brainparc`).

- hemis:

  character string vector, the hemispheres for each of the `vertices`.
  Allowed entries are `"lh"` and `"rh"`, for the left and right brain
  hemisphere, respectively. Must have same length as the 'vertices'
  vector (or exactly length 1, in which case we assume that this is the
  hemi for ALL query vertices).

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

## Value

a data.frame, the column names should be obvious.

## See also

[`coord_closest_regions`](https://dfsp-spirit.github.io/brainloc/reference/coord_closest_regions.md)
if you have a coordinate (on or near the surface) instead of a vertex.

## Examples

``` r
if (FALSE) { # \dontrun{
bp = brainparc_fs(get_subjects_dir(), "fsaverage", atlas="aparc");
vertex_closest_regions(bp, vertices=c(10, 20), hemis=c("lh", "rh"));
} # }
```
