# brainloc
Given a point on the human brain, identify standard space coordinates and find the closest brain regions according to a brain atlas.

This is currently intended to be used with FreeSurfer standard space templates (fsaverage, fsaverage6, ... fsaverage3).

## About

This is an R package that takes as input a vertex index of a FreeSurfer brain mesh in MNI305 space (typically fsaverage) and identifies the location in different coordinate systems and with respect to different brain regions. We use this to report the exact locations of clusters or other differences we find.

Supported coordinate systems include:

* MNI305 space RAS coordinates
* FreeSurfer 'MNI talairach' coordinates
* MNI152 space coordinates (using the linear transformation method, see [github.com/dfsp-spirit/regfusionr](https://github.com/dfsp-spirit/regfusionr) if you need a more accurate mapping).
* WIP: Talairach coordinates using Matthew Brett's transform from MNI152

The package also computes the closest brain regions and the distance to them, based on an atlas. You can use any atlas you like, the three default ones that come with FreeSurfer are:

* The Desikan-Killiany atlas (Desikan *et al.*, 2006. Neuroimage, 31(3):968-80)
* The Destrieux atlas (Destrieux *et al.*, 2010. Neuroimage, 53(1):1-15)
* The DKT40 altas from the [Mindboggle data set](https://mindboggle.info/data.html).

To use a different or custom atlas, just drop the respective annot files for the two hemispheres into the `/label/` directory of your template subject. See [the FreeSurfer documentation on Cortical Parcellations](https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation) for details on FreeSurfer brain atlases.



