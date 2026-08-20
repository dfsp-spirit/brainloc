# Download the talairach volume and labels from talairach.org.

Download the talairach volume `talairach.nii` and labels `labels.txt`
from talairach.org.

## Usage

``` r
download_talairach(accept_talairach_usage = TRUE)
```

## Value

Named list. The list has entries: "available": vector of strings. The
names of the files that are available in the local file cache. You can
access them using get_optional_data_file(). "missing": vector of
strings. The names of the files that this function was unable to
retrieve.

## Note

This function requires and internet connection. Files will only be
downloaded if they are not already available on the local computer.
