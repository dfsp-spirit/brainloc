# Create a hemilist of fs.annot instances from the given cluster overlay.

Create a
[`hemilist`](https://dfsp-spirit.github.io/brainloc/reference/hemilist.md)
of `fs.annot` instances from the given cluster overlay. An `fs.annot`
instance represents a brain surface parcellation, based on a brain
atlas. See
[`read.fs.annot`](https://rdrr.io/pkg/freesurferformats/man/read.fs.annot.html)
for details.

## Usage

``` r
clusteroverlay_to_annot(clusteroverlay, background_code = 0L, hemi = NULL)
```

## Arguments

- clusteroverlay:

  hemilist of integer vectors or a single integer vector of cluster
  overlay data: a vector that assigns each vertex to an integer class,
  all vertices belonging to a cluster share the same integer assignment.
  Non-cluster vertices are assigned the background class, see parameter
  'background_code'. One can also pass character strings, which will be
  interpreted as path to files that should be loaded with
  [`freesurferformats::read.fs.morph`](https://rdrr.io/pkg/freesurferformats/man/read.fs.morph.html)
  to get the clusteroverlay. One can also pass a `clusterinfo` instance,
  the 'overlay' field will be extracted and used in that case.

- background_code:

  scalar integer, the code in the overlayID data that should be
  interpreted as background or 'not part of any cluster'.

- hemi:

  character string, ignore unless clusteroverlay is a vector. In that
  case, it must be 'lh' or 'rh'. Used as a prefix when naming the
  clusters.

## Value

[`hemilist`](https://dfsp-spirit.github.io/brainloc/reference/hemilist.md)
of `freesurferformats::fs.annot` instances, each one representing a
brain surface parcellation for one hemisphere. Each label in the
parcellation (with the exception of the unknown/background region)
corresponds to a cluster.
