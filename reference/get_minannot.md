# Extract minannot from fs.annot instance.

Extract minannot from fs.annot instance.

## Usage

``` r
get_minannot(fs.annot)
```

## Arguments

- fs.annot:

  a `freesurferformats::fs.annot` instance, or a `minannot` instance,
  which will be returned as-is. Alternatively hemilist of the latter.

## Value

a `minannot` instance, i.e., only the 'label_names' field of the
`fs.annot` instance.
