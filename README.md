# brainloc
Given a vertex (point) or vertex set (cluster) on the human brain, identify standard space coordinates and find the closest and/or overlapping brain regions according to a brain atlas.

This is currently intended to be used with [FreeSurfer](https://freesurfer.net/) standard space templates (fsaverage, fsaverage6, ... fsaverage3).


![Fig1](./web/brainloc.png?raw=true "Brainloc.")
**Fig. 1** *Left: Vertex 145029 on the left fsaverage surface (at pink marker). Screenshot from the FreeView application that comes with [FreeSurfer](https://freesurfer.net). Right: Location of MNI coordinate 39 -30  65, the result of mapping fsaverage vertex 145029 to MNI152 space. Screenshot from the [MNI - Talairach Tool](https://bioimagesuiteweb.github.io/bisweb-manual/tools/mni2tal.html).* 


## Features

This is an R package that takes as input a vertex index of a FreeSurfer brain mesh in MNI305 space (typically fsaverage) and identifies the location in different coordinate systems and with respect to different brain regions. We use this to report the exact locations of clusters or other differences we find. It can also be used to find the vertex in a brain mesh that is closest to a given coordinate.

### Full list of features and functions implementing them

- Given a coordinate in MNI305 space:
  - `coord_closest_vertex()`: find the closest vertex in a mesh that is MNI305 space (e.g., fsaverage, fsaverage6, etc).
  - Note: This actually works for any brain mesh and a coordinate in its surface space, including native space meshes and any surface (white, pial, ...).
- Given a vertex or coordinate on a mesh in MNI305 space:
  - `coord_MNI305_info()`: find its MNI305 coordinate (trivial)
  - `coord_MNI305_info()`: find its MNI152 coordinate using one of the following two methods:
    - with the FreeSurfer 4x4 matrix or
    - with the more accurate [regfusionr](https://github.com/dfsp-spirit/regfusionr) method by [Wu *et al.*](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6239990/) (requires optional [regfusionr](https://github.com/dfsp-spirit/regfusionr) package).
  - `coord_MNI305_info()`: find its Talairach coordinates (using [Matthew Brett's non-linear transform from MNI152](https://brainmap.org/training/BrettTransform.html))
  - `vertex_closest_regions()`: find the region the vertex is assigned to in a brain atlas parcellation like the Desikan atlas that comes with FreeSurfer (trivial)
  - `vertex_closest_regions()`: find the distances to all other atlas regions, with different distance methods (Euclidean, geodesic along the mesh) and different linkages (defining the reference point when measuring the distance to a region, e.g., closest vertex in region, or center vertex)
  - `coord_closest_regions()`: do the same for a coordinate instead of a vertex index.
- Given a cluster as a set of vertices on any brain surface mesh:
  - `clusterinfo()`: read such cluster information from a statistical map (e.g., t-value map for all mesh vertices) and an overlay map assigning a cluster identifier to each vertex.
  - `clusterinfo_from_thresholded_overlay()`: read such cluster information only from a thresholded statistical map, using BFS on the mesh to identify the clusters.
  - `cluster_location_details()`: find the extremum value and vertex of each cluster
  - `cluster_peaks()`: find all peaks of each cluster
  - `cluster_region_overlap()`: given an additional brain parcellation, find all regions the cluster overlaps with and compute the percentage overlap for both the cluster and the regions.
- Given a coordinate in Talairach space:
  - `get_talairach_label()`: retrieve the Talairach volume label for the point (5 level hierarchy, e.g., 'Right Cerebrum, Temporal Lobe, Sub-Gyral, Gray Matter, Brodmann area 20'). Please cite the two [Talairach References papers listed on talairach.org](http://www.talairach.org/) when using this functionality.
  
 
## Documentation

* A detailed vignette with explanations and examples for the functions of the package is included, run `browseVignettes("brainloc")` to see the vignette. If the last build succeeded, you may be able to [read the vignette online here on r-universe](https://dfsp-spirit.r-universe.dev/ui#view:brainloc/brainloc.html).
* Help for a specific function can be accessed in the usual R manner: `?<function>`, where you replace `<function>` with a function name. Like this: `?vertex_closest_regions`.
* To see all functions which are part of the package API, run: `help(package="brainloc")`.
* Run `example(<function>)` to see a live demo that uses the function `<function>`. Like this: `example(vertex_closest_regions)`.
* The [unit tests](./tests/testthat/) that come with this package are essentially a list of examples that illustrate how to use the functions.


## Demonstrations

The source code that was used to generate these images is [available in this unit test file](tests/testthat/test-visualization_with_fsbrain.R).


![Fig2](./web/brainloc_coordinate_closest_vertex.png?raw=true "Brainloc: closest vertex to a query coordinate.")

**Fig. 2** *An arbitrary query coordinate (red sphere) and the closest mesh vertex (yellow sphere).*

![Fig3](./web/brainloc_parcellation_region_center_neighborhoods.png?raw=true "Brainloc: brain atlas regions with neighborhood relationships.")

**Fig. 3** *Brain atlas regions with neighborhood relationship. Each sphere representes the centroid of a brain atlas region from the Desikan atlas. Each edge represents a spatial neighborhood relationship between the respective pair of brain regions. Region data from the fsaverage subject.*



## Installation

It's still rather early, but if you want to try the current version, run the following commands from an `R` session:

```R
options(repos = c(
    dfspspirit = 'https://dfsp-spirit.r-universe.dev',
    CRAN = 'https://cloud.r-project.org'))

install.packages('brainloc')
```
