# Given vertex indices defining a cluster, find all atlas regions the cluster overlaps with.

Given vertex indices defining a cluster, find all atlas regions the
cluster overlaps with.

## Usage

``` r
cluster_overlapping_regions(annot_min, cluster_vertices)
```

## Arguments

- annot_min:

  a full `fs.annot` instance, or a minimal annot, i.e., only the
  'label_names' field of the annot. Only a single one, not a
  [`hemilist`](https://dfsp-spirit.github.io/brainloc/reference/hemilist.md).

- cluster_vertices:

  integer vector, the vertices defining the cluster (technically they do
  not need to form a cluster or be connected). Must not be empty.

## Value

data.frame with columns 'region': the region name,
'num_shared_vertices': the number of cluster vertices which are in the
region, and 'percent_shared_vertices': the percent of cluster vertices
which are in the region, 'cluster_percent_of_region': how much of the
region area is filled by the cluster vertices, in percent. The
data.frame rows are ordered descending by 'percent_shared_vertices'.
