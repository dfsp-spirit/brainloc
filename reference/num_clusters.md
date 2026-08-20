# Compute number of clusters from cluster annot or clusterinfo instance.

Compute number of clusters from cluster annot or clusterinfo instance.

## Usage

``` r
num_clusters(cluster_annots)
```

## Arguments

- cluster_annots:

  [`hemilist`](https://dfsp-spirit.github.io/brainloc/reference/hemilist.md)
  of cluster annots, see
  [`clusteroverlay_to_annot`](https://dfsp-spirit.github.io/brainloc/reference/clusteroverlay_to_annot.md).

## Value

named list with keys 'lh', 'rh' and 'total', each holding a scalar
integer. The cluster counts.
