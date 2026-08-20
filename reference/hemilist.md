# Create hemilist.

A hemilist is a named list which has at least one of the following
fields: `lh` and `rh`. The fields can store arbitrary values.

## Usage

``` r
hemilist(lh_data = NULL, rh_data = NULL)
```

## Arguments

- lh_data:

  any R object, the lh value.

- rh_data:

  any R object, the rh value.

## Value

a hemilist

## Note

This function is trivial and mainly exists to document what a hemilist
is (so this explanation can be linked to from the documentation of other
functions which use the concept).

## See also

Other hemilist functions:
[`hemilist.wrap()`](https://dfsp-spirit.github.io/brainloc/reference/hemilist.wrap.md),
[`is.hemilist()`](https://dfsp-spirit.github.io/brainloc/reference/is.hemilist.md)
