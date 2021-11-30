# brainloc
Given a vertex (point) or vertex set (cluster) on the human brain, identify standard space coordinates and find the closest brain regions according to a brain atlas.

This is currently intended to be used with FreeSurfer standard space templates (fsaverage, fsaverage6, ... fsaverage3).


![Fig1](./web/brainloc.png?raw=true "Brainloc.")

## Features

This is an R package that takes as input a vertex index of a FreeSurfer brain mesh in MNI305 space (typically fsaverage) and identifies the location in different coordinate systems and with respect to different brain regions. We use this to report the exact locations of clusters or other differences we find. It can also be used to find the vertex in a brain mesh that is closest to a given coordinate.


1) Given a coordinate in MNI305 space:
 - find the closest vertex in a mesh that is MNI305 space (e.g., fsaverage, fsaverage6, etc).
2) Given a vertex on a mesh in MNI305 space:
 - find its MNI305 coordinate (trivial)
 - find its MNI152 coordinate
   - using the FreeSurfer 4x4 matrix or
   - using the more accurate [regfusionr](https://github.com/dfsp-spirit/regfusionr) method
 - find its Talairach coordinates (using Matthew Brett's non-linear transform from MNI152)
 - find the region the vertex is assigned to in a brain atlas parcellation like the Desikan atlas that comes with FreeSurfer (trivial)
 - find the distances to all other atlas regions, with different distance methods (Euclidean, geodesic along the mesh) and different linkages (defining the reference point when measuring the distance to a region, e.g., closest vertex in region, or center vertex)
3) Given a cluster as a set of vertices on any brain surface mesh:
 - read such cluster information from a statistical map (e.g., t-value map for all mesh vertices) and an overlay map assigning a cluster identifier to each vertex.
 - read such cluster information only from a thresholded statistical map, using BFS on the mesh to identify the clusters.
 - find the extremum value of each cluster
 - find all peaks of each cluster
 - given a brain parcellation, find all regions the cluster overlaps with and compute the percentage overlap for both the cluster and the regions
 


## Installation

Not yet, this is work-in-progress. Come back another day.

