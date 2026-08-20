# Get details on cluster locations.

Get details on cluster locations.

## Usage

``` r
cluster_location_details(
  clusterinfo,
  silent = getOption("brainloc.silent", default = FALSE),
  ...
)
```

## Arguments

- clusterinfo:

  a clusterinfo instance, see the `clusterinfo` function on how to get
  one. Must contain a valid `brainparc` in field 'brainparc'.

- silent:

  logical, whether to suppress console messages.

- ...:

  passed on to
  [`cluster_extrema`](https://dfsp-spirit.github.io/brainloc/reference/cluster_extrema.md).

## Value

a data.frame with cluster location details, including MNI152 and
Talairach coordinates, Talairach Daemon labels (5 hierarchy levels +
full label), and the atlas region name from the first atlas in the
brainparc. The column names should be self-explanatory.
