# Find subjects_dir or stop.

Try a number of ways to find the subjects dir, using environment
variables and knowledge on common installation paths of FreeSurfer.

## Usage

``` r
get_subjects_dir(
  mustWork = TRUE,
  allow_download = FALSE,
  accept_freesurfer_license = FALSE
)
```

## Arguments

- mustWork:

  whether to stop if no subjects_dir can be found.

- allow_download:

  logical, whether to allow downloading the data in case it is not found
  locally.

- accept_freesurfer_license:

  logical, whether you want to also download fsaverage and fsaverage3,
  and accept the FreeSurfer license for fsaverage and fsaverage3,
  available at
  https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense.
  Defaults to FALSE. If FALSE, nothing will be downloaded.

## Value

character string, the subjects directory. If mustWork is FALSE, it will
return NULL if no directory was found. If mustWork is TRUE and nothing
is found, it stops.

## See also

Other FreeSurfer helper functions:
[`find.freesurferhome()`](https://dfsp-spirit.github.io/brainloc/reference/find.freesurferhome.md),
[`fs_home()`](https://dfsp-spirit.github.io/brainloc/reference/fs_home.md),
[`has_fs()`](https://dfsp-spirit.github.io/brainloc/reference/has_fs.md)
