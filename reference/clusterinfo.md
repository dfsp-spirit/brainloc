# Create clusterinfo data structure from files or pre-loaded data.

Create clusterinfo data structure from files or pre-loaded data.

## Usage

``` r
clusterinfo(
  lh_overlay,
  rh_overlay,
  lh_statmap,
  rh_statmap,
  template_subject = "fsaverage",
  subjects_dir = file.path(getOption("brainloc.fs_home", default =
    Sys.getenv("FREESURFER_HOME")), "subjects")
)
```

## Arguments

- lh_overlay:

  integer vector, the cluster overlay data: one integer per vertex of
  the left hemisphere. This assigns each vertex to a cluster, and all
  vertices of one cluster have the same number. Background is typically
  1L. If a character string, the parameter will be interpreted as a file
  path and loaded with
  [`freesurferformats::read.fs.morph`](https://rdrr.io/pkg/freesurferformats/man/read.fs.morph.html).

- rh_overlay:

  integer vector, just like `lh_overlay`, but for the right hemisphere.

- lh_statmap:

  double vector, the stats map. Typically a t-value map. Must have one
  value per vertex. If a character string, the parameter will be
  interpreted as file path and loaded with
  [`freesurferformats::read.fs.morph`](https://rdrr.io/pkg/freesurferformats/man/read.fs.morph.html).

- rh_statmap:

  double vector, just like `lh_statmap`, but for the right hemisphere.

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

named list with entries 'overlay', 'statmap', and 'metadata': a
clusterinfo data structure. Each of the 'overlay' and 'statmap' keys
holds a
[`hemilist`](https://dfsp-spirit.github.io/brainloc/reference/hemilist.md)
of numerical vectors.

## See also

[`clusterinfo_from_thresholded_overlay`](https://dfsp-spirit.github.io/brainloc/reference/clusterinfo_from_thresholded_overlay.md)
If you do have a thresholded t-map instead of one t-map and one cluster
overlay map.
