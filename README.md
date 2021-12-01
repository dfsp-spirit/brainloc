# brainloc
Given a vertex (point) or vertex set (cluster) on the human brain, identify standard space coordinates and find the closest and/or overlapping brain regions according to a brain atlas.

This is currently intended to be used with FreeSurfer standard space templates (fsaverage, fsaverage6, ... fsaverage3).


![Fig1](./web/brainloc.png?raw=true "Brainloc.")

## Features

This is an R package that takes as input a vertex index of a FreeSurfer brain mesh in MNI305 space (typically fsaverage) and identifies the location in different coordinate systems and with respect to different brain regions. We use this to report the exact locations of clusters or other differences we find. It can also be used to find the vertex in a brain mesh that is closest to a given coordinate.

### Full list of features and functions implementing them

- Given a coordinate in MNI305 space:
  - `coord_closest_vertex()`: find the closest vertex in a mesh that is MNI305 space (e.g., fsaverage, fsaverage6, etc).
  - Note: This actually works for any brain mesh and a coordinate in its surface space, including native space meshes and any surface (white, pial, ...).
- Given a vertex or coordinate on a mesh in MNI305 space:
  - `coord_MNI305_info()`: find its MNI305 coordinate (trivial): see function `coord_MNI305_info()`
  - `coord_MNI305_info()`: find its MNI152 coordinate using one of the following two methods:
    - with the FreeSurfer 4x4 matrix or
    - with the more accurate [regfusionr](https://github.com/dfsp-spirit/regfusionr) method by Wu *et al.* (requires optional regfusionr package).
  - `coord_MNI305_info()`: find its Talairach coordinates (using Matthew Brett's non-linear transform from MNI152) `coord_MNI305_info()`
  - `vertex_closest_regions()`: find the region the vertex is assigned to in a brain atlas parcellation like the Desikan atlas that comes with FreeSurfer (trivial)
  - `vertex_closest_regions()`: find the distances to all other atlas regions, with different distance methods (Euclidean, geodesic along the mesh) and different linkages (defining the reference point when measuring the distance to a region, e.g., closest vertex in region, or center vertex)
  - `coord_closest_regions()`: do the same for a coordinate instead of a vertex index.
- Given a cluster as a set of vertices on any brain surface mesh:
  - `clusterinfo()`: read such cluster information from a statistical map (e.g., t-value map for all mesh vertices) and an overlay map assigning a cluster identifier to each vertex.
  - `clusterinfo_from_thresholded_overlay()`: read such cluster information only from a thresholded statistical map, using BFS on the mesh to identify the clusters.
  - `cluster_location_details()`: find the extremum value and vertex of each cluster
  - `cluster_peaks()`: find all peaks of each cluster
  - `cluster_region_overlap()`: given an additional brain parcellation, find all regions the cluster overlaps with and compute the percentage overlap for both the cluster and the regions.
 


## Installation

It's still a bit early, but if you want to try the current version, run the following commands from an `R` session:

```R
options(repos = c(
    dfspspirit = 'https://dfsp-spirit.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))

install.packages('brainloc')
```
