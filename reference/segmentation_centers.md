# Find center of mass of regions in a brain volume segmentation.

Find center of mass of regions in a brain volume segmentation.

## Usage

``` r
segmentation_centers(volatlas, named_regions = NULL, vox2ras = diag(4))
```

## Arguments

- volatlas:

  3d integer array, the brain volume segmentation. Can be an `fs.volume`
  instance. Typically something like `'subject/mri/aseg.mgz'`.

- named_regions:

  the regions to consider, a named list where the keys are character
  strings: the region names from the colorLUT, and the values are
  integers: the region codes from the volume. See
  [`name_regions`](https://dfsp-spirit.github.io/brainloc/reference/name_regions.md)
  to get one. Can be `NULL`, in which case all regions in the volatlas
  will be used, and the region names will appear as 'region_x', where
  'x' is the integer region code in the volatlas.

- vox2ras:

  4x4 numerical matrix, the voxel to RAS coordinate transformation
  matrix for the volume. If NULL and volatlas is an fs.volume instance,
  the function will try to extract it from the volume header. The
  default is the identity matrix, which will return voxel indices
  instead of RAS coordinates.

## Value

data.frame with region names and the center of mass x, y and z
coordinates for the regions.

## Examples

``` r
if (FALSE) { # \dontrun{
if(has_fs()) {
fsh = fs_home();
lutfile = file.path(fsh, "FreeSurferColorLUT.txt");
segfile = file.path(fsh, 'subjects', 'fsaverage', 'mri', 'aseg.mgz');
named_regions = name_regions(segfile, lutfile);
named_regions$Unknown = NULL; # Ignore an atlas region.
sc = segmentation_centers(segfile, named_regions);
regions_distance_matrix = dist(sc);
}
} # }
```
