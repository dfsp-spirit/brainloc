# Extract the vertex coordinates from the brainparc surface.

Extract the vertex coordinates from the brainparc surface.

## Usage

``` r
get_surface_coords(brainparc, vertices, hemis)
```

## Arguments

- brainparc:

  a `brainparc` instance, see
  [`brainparc_fs`](https://dfsp-spirit.github.io/brainloc/reference/brainparc_fs.md)
  to get one.

- vertices:

  integer vector, the query vertex indices.

- hemis:

  vector of character string, the hemispheres for the query vertices.
  Length must match, length of 1 will be recycled to all vertices.

## Value

numerical nx3 matrix, the coordindates for the x query vertices.

## See also

Other brainparc accessors:
[`get_surface()`](https://dfsp-spirit.github.io/brainloc/reference/get_surface.md)
