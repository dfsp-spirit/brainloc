# Compute overlayID map from thresholded statmap (clustermap).

This uses breadth-first search on the graph of the mesh to identify all
connected cluster values and identify the separate clusters. Different
clusters must not touch for this to work (if they do, there is no
reliable way to obtain the clusters from the thresholded map). Clusters
can have holes.

## Usage

``` r
clusteroverlay_from_threshmap(threshmap, value_thresholded, surface)
```

## Arguments

- threshmap:

  double vector, the stats map. Typically a thresholded t-value map
  (cluster map). Must have one value per vertex. The value assigned to
  vertices that have been removed by the thresholding can be set with
  parameter 'value_thresholded'.

- value_thresholded:

  scalar double, the data value of thresholded vertices in the
  `lh_threshmap` and `rh_threshmap`.

- surface:

  a single `fs.surface` instance, used for neighborhood computation.

## Value

integer vector, the clusterID overlay.

## Note

This visits every node and edge in the mesh/graph exactly once.
