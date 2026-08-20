# Get color string for a region code in a color lookup table.

Get color string for a region code in a color lookup table.

## Usage

``` r
region_color(colorlut, struct_idx)
```

## Arguments

- colorlut:

  a color lookup table data.frame, as returned by
  [`read.fs.colortable`](https://rdrr.io/pkg/freesurferformats/man/read.fs.colortable.html).
  Can be `NULL` if you do not have any, in that case the region names
  will just be 'region_x', where 'x' is the integer region code in the
  volatlas.

- struct_idx:

  integer vector, the region / structure codes you are interested in.

## Value

vector of character strings, the hex color names (like '#FF0000' for
red).
