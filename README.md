# brainloc
Given a point on the human brain, identify standard space coordinates and find the closest brain regions according to a brain atlas.

## About

This is an R package that takes as input a vertex index of a FreeSurfer brain mesh in MNI305 space (typically fsaverage) and identifies the location in different coordinate systems, including:

* MNI305 space RAS coordinates
* FreeSurfer 'MNI talairach' coordinates
* MNI152 space coordinates.

The package also computes the closest brain regions, based on an atlas.
