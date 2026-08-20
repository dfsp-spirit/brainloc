# Get the clusters from a clusterinfo instance.

Get the clusters from a clusterinfo instance.

## Usage

``` r
get_clusters(clusterinfo, hemi = "both")
```

## Arguments

- clusterinfo:

  a `clusterinfo` instance.

- hemi:

  character string, one of 'lh', 'rh' or 'both'. The hemisphere(s) for
  which to return the clusters.

## Value

named list, the keys are the cluster names, and the values are integer
vectors defining the member vertices.
