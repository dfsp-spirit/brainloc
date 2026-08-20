# Load a surface for a subject.

Load a brain surface mesh for a subject.

## Usage

``` r
subject.surface(
  subjects_dir,
  subject_id,
  surface = "white",
  hemi = "both",
  force_hemilist = FALSE
)
```

## Arguments

- subjects_dir:

  string. The FreeSurfer SUBJECTS_DIR, i.e., a directory containing the
  data for all your subjects, each in a subdir named after the subject
  identifier.

- subject_id:

  string. The subject identifier

- surface:

  string. The surface name. E.g., "white", or "pial". Used to construct
  the name of the surface file to be loaded.

- hemi:

  string, one of 'lh', 'rh', or 'both'. The hemisphere name. Used to
  construct the names of the surface file to be loaded. For 'both', see
  the information on the return value.

- force_hemilist:

  logical, whether to return a
  [`hemilist`](https://dfsp-spirit.github.io/brainloc/reference/hemilist.md)
  even if the 'hemi' parameter is not set to 'both'

## Value

the \`fs.surface\` instance, as returned by
[`read.fs.surface`](https://rdrr.io/pkg/freesurferformats/man/read.fs.surface.html).
If parameter \`hemi\` is set to \`both\`, a named list with entries
\`lh\` and \`rh\` is returned, and the values of are the respective
surfaces. The mesh data structure used in \`fs.surface\` is a \*face
index set\*.
