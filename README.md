# brainloc
Given a point on the human brain, identify standard space coordinates and find the closest brain regions according to a brain atlas.

This is currently intended to be used with FreeSurfer standard space templates (fsaverage, fsaverage6, ... fsaverage3).

## About

This is an R package that takes as input a vertex index of a FreeSurfer brain mesh in MNI305 space (typically fsaverage) and identifies the location in different coordinate systems and with respect to different brain regions. We use this to report the exact locations of clusters or other differences we find.

Supported coordinate systems include:

* MNI305 space RAS coordinates (simply the coordinate of the input vertex index)
* MNI152 space coordinates (using the linear transformation method with the FreeSurfer matríx ([section 8b here](https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems)), see [github.com/dfsp-spirit/regfusionr](https://github.com/dfsp-spirit/regfusionr) if you need a more accurate mapping).
* Talairach coordinates using [Matthew Brett's transform](http://brainmap.org/training/BrettTransform.html) from MNI152

The package also computes the closest brain regions and the distance to them, based on an atlas. You can use any atlas you like, the three default ones that come with FreeSurfer are:

* The Desikan-Killiany atlas (Desikan *et al.*, 2006. Neuroimage, 31(3):968-80)
* The Destrieux atlas (Destrieux *et al.*, 2010. Neuroimage, 53(1):1-15)
* The DKT40 altas from the [Mindboggle data set](https://mindboggle.info/data.html).

To use a different or custom atlas, just drop the respective annot files for the two hemispheres into the `/label/` directory of your template subject. See [the FreeSurfer documentation on Cortical Parcellations](https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation) for details on FreeSurfer brain atlases.


### Atlas distances

The following methods are available to compute the distance of a point on a brain surface to a brain atlas region:

* Euclidean distance: point to mean value of region coordinates
* Euclidean distance: point to closest vertex of region

The following methods are work in progress:

* Geodesic distance to closest vertex of region

