# Compute adjacency matrix of surface parcellation (annot) regions.

Given a brain surface and a parcellation into regions, find out which
regions are adjacent to each other. A region `i` is adjacent to another
region `j` if an edge connects any vertex of `i` with any vertex of `i`.

## Usage

``` r
annot_neighbors(annot_min, surface, empty_rename = "_")
```

## Arguments

- annot_min:

  a full `fs.annot` instance, or a minimal annot, i.e., only the
  'label_names' field of the annot. Only a single one, not a
  [`hemilist`](https://dfsp-spirit.github.io/brainloc/reference/hemilist.md).

- surface:

  `fs.surface` instance, the brain mesh (for a single hemisphere).

- empty_rename:

  optional character string, the name to use for unnamed regions (which
  have an empty string as region name).

## Value

named integer matrix of regions, expressing whether they are direct
neighbors (value `1L`) or not (value `0L`).

## Note

Regions from different hemispheres cannot be direct neighbors.

This works on a single hemisphere. We should create a wrapper that works
on brainparc instances.
