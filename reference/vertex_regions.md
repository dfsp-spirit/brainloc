# Find regions for vertices in all brainparc atlases.

Find regions for vertices in all brainparc atlases.

## Usage

``` r
vertex_regions(brainparc, vertices, hemis)
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

## Value

a data.frame, the column names should be obvious.

## See also

[`vertex_closest_regions`](https://dfsp-spirit.github.io/brainloc/reference/vertex_closest_regions.md)
to compute distance to all regions.

## Examples

``` r
if (FALSE) { # \dontrun{
bp = brainparc_fs(get_subjects_dir(), "fsaverage", atlas="aparc");
vertex_region(bp, vertices=c(10, 20), hemis=c("lh", "rh"));
} # }
```
