# Download FreeSurfer template and demo data if needed and return its path (subjects_dir).

This is a wrapper around
[`download_fsaverage()`](https://dfsp-spirit.github.io/brainloc/reference/download_fsaverage.md)
and `download_fsaverage3())`. It will download the data from the
internet unless it already exists locally.

## Usage

``` r
sjd_demo(accept_freesurfer_license = FALSE)
```

## Arguments

- accept_freesurfer_license:

  logical, whether you want to also download fsaverage and fsaverage3,
  and accept the FreeSurfer license for fsaverage and fsaverage3,
  available at
  https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense.
  Defaults to FALSE. If FALSE, nothing will be downloaded.

## Value

character string, the path to the 'subjects_dir' directory within the
downloaded template data directory.

## Note

This function will stop if the data cannot be accessed, i.e., the
'subjects_dir' does not exist after trying to download the data.
