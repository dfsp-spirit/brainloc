# brainloc
Given a point on the human brain, identify standard space coordinates and find the closest brain regions according to a brain atlas.

This is currently intended to be used with FreeSurfer standard space template (fsaverage, fsaverage6, ... fsaverage3).

## About

This is an R package that takes as input a vertex index of a FreeSurfer brain mesh in MNI305 space (typically fsaverage) and identifies the location in different coordinate systems, including:

* MNI305 space RAS coordinates
* FreeSurfer 'MNI talairach' coordinates
* MNI152 space coordinates.

The package also computes the closest brain regions, based on an fsaverage atlas. You can use any atlas you like, the three default ones that come with FreeSurfer for fsaverage are:

* The Desikan-Killiany atlas (Desikan *et al.*, 2006. Neuroimage, 31(3):968-80)
* The Destrieux atlas (Destrieux *et al.*, 2010. Neuroimage, 53(1):1–15)
* DKT40 altas from the [Mindboggle data set](https://mindboggle.info/data.html)

See [the FreeSurfer documentation on Cortical Parcellations](https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation) for details.
