# Return FreeSurfer path.

Return FreeSurfer path.

## Usage

``` r
fs_home()
```

## Value

the FreeSurfer path, typically what the environment variable
\`FREESURFER_HOME\` points to.

## Note

This function will stop (i.e., raise an error) if the directory cannot
be found. It calls
[`find.freesurferhome`](https://dfsp-spirit.github.io/brainloc/reference/find.freesurferhome.md)
internally, see there for more details.

## See also

Other FreeSurfer helper functions:
[`find.freesurferhome()`](https://dfsp-spirit.github.io/brainloc/reference/find.freesurferhome.md),
[`get_subjects_dir()`](https://dfsp-spirit.github.io/brainloc/reference/get_subjects_dir.md),
[`has_fs()`](https://dfsp-spirit.github.io/brainloc/reference/has_fs.md)
