# Transform FreeSurfer surface space coordinates into MNI Talairach space using talairach.xfm file.

This uses the affine transform in
`<subject>/mri/transforms/talairach.xfm`, be sure to know about the
limitations of this approach. See the FreeSurfer documentation on
talairach for details. The coordinates produced by this approach are
referred to as 'MNI Talairach' coordinates in the FreeSurfer
documentation, and are not identical to Talairach space.

## Usage

``` r
coord_fssurface_to_mnitalairach(subjects_dir, subject_id, surface_coords)
```

## Arguments

- subjects_dir:

  character string, the path to the recon-all SUBJECTS_DIR containing
  the data of all subjects of your study.

- subject_id:

  character string, the subject identifier (i.e., the sub directory name
  under the `subjects_dir`). This can be a template like fsaverage, but
  it also works for other subjects.

- surface_coords:

  nx3 numerical matrix of surface coordinates (vertex positions in
  FreeSurfer surface space).

## Value

nx3 numerical matrix of MNI Talairach space coordinates.

## Examples

``` r
if (FALSE) { # \dontrun{
# Will lead to identity transform for fsaverage, so no change in coords.
sjd = fsbrain::fsaverage.path();
sj = "fsaverage";
surf = fsbrain::subject.surface(sjd, sj, hemi="lh", surface="orig");
tal = coord_fssurface_to_mnitalairach(sjd, sj, surf$vertices);
} # }
```
