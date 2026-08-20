# Compute cluster peaks.

A peak is a position in a cluster (a vertex) which is assigned the
maximum statsmap value over its neighborhood. The neighborhood
definition used by this function is the 1-ring neighborhood on the mesh.

## Usage

``` r
cluster_peaks(
  clusterinfo,
  type = "extreme",
  silent = getOption("brainloc.silent", default = FALSE),
  ...
)
```

## Arguments

- clusterinfo:

  a clusterinfo instance, see the `clusterinfo` function to get one.

- type:

  character string, one of 'extreme', 'min' or 'max'. Which cluster
  value to report. The default of 'extreme' reports the absolutely
  larger one of the min and the max value and is for the typical use
  case of `t`-value maps.

- silent:

  whether to suppress printing messages to stdout.

- ...:

  passed on to
  [`clusteroverlay_to_annot`](https://dfsp-spirit.github.io/brainloc/reference/clusteroverlay_to_annot.md).

## Value

a `data.frame` with cluster peak information. The column names should be
self-explanatory.
