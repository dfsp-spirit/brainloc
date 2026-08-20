# Compute peaks of a single cluster.

Compute peaks of a single cluster.

## Usage

``` r
single_cluster_peaks(cluster_vertices, statmap, surface, type = "extreme")
```

## Arguments

- cluster_vertices:

  integer vector, the vertex indices of the cluster.

- statmap:

  double vector, the full stat map for the surface. Only the values of
  the cluster_vertices are used.

- surface:

  a single `fs.surface` instance. Used to compute the neighborhood of
  the cluster vertices.

- type:

  character string, one of 'extreme', 'min' or 'max'. Which cluster
  value to report. The default of 'extreme' reports the absolutely
  larger one of the min and the max value and is for the typical use
  case of `t`-value maps.

## Value

named list with entries 'vertex' and 'value', they contain an integer
vector and a double vector, respectively.

## Note

Called by `cluster_peaks`, use that instead.
