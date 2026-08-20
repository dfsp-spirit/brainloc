# Compute adjacency of surface parcellation (annot) regions.

Given a brain surface and a parcellation into regions, find out which
regions are adjacent to each other. A region `i` is adjacent to another
region `j` if an edge connects any vertex of `i` with any vertex of `i`.
Works with brainparcellations.

## Usage

``` r
brainparc_neighbors(bp)
```

## Arguments

- bp:

  a brainparc instance

## Value

hemilist of named integer matrix of regions, expressing whether they are
direct neighbors (value `1L`) or not (value `0L`).
