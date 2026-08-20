# Create clusterinfo instance from thresholded per-vertex data overlays.

Create clusterinfo instance from thresholded per-vertex data overlays.

## Usage

``` r
clusterinfo_from_thresholded_overlay(
  lh_threshmap,
  rh_threshmap,
  value_thresholded = 0,
  template_subject = "fsaverage",
  subjects_dir = file.path(getOption("brainloc.fs_home", default =
    Sys.getenv("FREESURFER_HOME")), "subjects")
)
```

## Arguments

- lh_threshmap:

  double vector, the stats map. Typically a thresholded t-value map.
  Must have one value per vertex. The value assigned to vertices that
  have been removed by the thresholding can be set with parameter
  'value_thresholded'. If a character string, the parameter will be
  interpreted as file path and loaded with
  [`freesurferformats::read.fs.morph`](https://rdrr.io/pkg/freesurferformats/man/read.fs.morph.html).

- rh_threshmap:

  like `lh_threshmap`, but for the right hemisphere.

- value_thresholded:

  scalar double, the data value of thresholded vertices in the
  `lh_threshmap` and `rh_threshmap`.

- template_subject:

  character string, the template subject name. Typically 'fsaverage' or
  'fsaverage6'. Must be a subject in MNI305 space (i.e., the standard
  FreeSurfer template space). \*\*Important:\*\* The overlay and statmap
  data passed to this function must be in the same space as the
  template_subject (usually MNI305 / fsaverage space). Passing data from
  a different space (e.g., native subject space) will produce silently
  incorrect coordinate transforms and atlas region assignments.

- subjects_dir:

  character string, file system path to a directory containing the
  recon-all data for the template_subject. Used to load surfaces and
  annotations to identify cluster coordinates and atlas regions.

## Value

a `clusterinfo` instance.

## Note

This can only work if there are no clusters which touch each other. I
think this can never happen, but we should double-check.

## See also

[`clusterinfo`](https://dfsp-spirit.github.io/brainloc/reference/clusterinfo.md)
If you do have a separate t-map and a cluster overlay map instead of
only a thresholded t-map.
