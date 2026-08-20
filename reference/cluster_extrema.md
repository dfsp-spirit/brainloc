# Print information on cluster extreme values.

A cluster extreme value is the most extreme value of a cluster. One can
choose whether the minimum, maximum, or absolute maximum is requested.
Each cluster has exactly one extreme value (though it may occur at
several vertices in rare cases).

## Usage

``` r
cluster_extrema(
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

a `data.frame` with cluster extrema information. The column names should
be self-explanatory.

## Note

If the extreme value occurs at several vertices of a cluster, which of
these vertices will be reported is undefined.
