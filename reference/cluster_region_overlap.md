# Compute overlapping atlas regions for clusters.

For each cluster in the clusterinfo, compute 1) which atlas regions it
overlaps with, and to what percentage the cluster is contained in each
region. And, 2), how much of each of the overlapping regions are covered
by the cluster. The atlases or surface parcellations to compare the
clusters to are take from the `clusterinfo$brainparc` field.

## Usage

``` r
cluster_region_overlap(
  clusterinfo,
  silent = getOption("brainloc.silent", default = FALSE)
)
```

## Arguments

- clusterinfo:

  a clusterinfo instance, see the `clusterinfo` function to get one.

- silent:

  whether to suppress printing messages to stdout.

## Value

a data.frame with cluster overlap information.
