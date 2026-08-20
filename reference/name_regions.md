# Get names of regions in labeled volume from color lookup table (colorLUT).

Given a brain volume with integer labels as voxel values and a suitable
colorLUT, extract the region names as a named list.

## Usage

``` r
name_regions(
  volatlas,
  colorlut,
  ignore_not_in_lut = FALSE,
  warn_not_in_lut = FALSE
)
```

## Arguments

- volatlas:

  3d integer array, the brain volume segmentation. Can be an `fs.volume`
  instance. Typically something like `'subject/mri/aseg.mgz'`.

- colorlut:

  a color lookup table data.frame, as returned by
  [`read.fs.colortable`](https://rdrr.io/pkg/freesurferformats/man/read.fs.colortable.html).
  Can be `NULL` if you do not have any, in that case the region names
  will just be 'region_x', where 'x' is the integer region code in the
  volatlas.

- ignore_not_in_lut:

  logical, whether to ignore regions that are not found in the colorLUT.
  If TRUE, the regions will not appear in the named output list. If
  FALSE, they will appear as 'region_x', where 'x' is the integer region
  code in the volatlas.

- warn_not_in_lut:

  logical, whether to print a warning if some regions of the volume do
  not occur in the colorLUT.

## Value

named list, the keys are character strings: the region names from the
colorLUT, and the values are integers: the region codes from the volume.
