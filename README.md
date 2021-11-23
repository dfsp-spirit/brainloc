# brainloc
Given a point on the human brain, identify standard space coordinates and find the closest brain regions according to a brain atlas.

This is currently intended to be used with FreeSurfer standard space templates (fsaverage, fsaverage6, ... fsaverage3).

## About

This is an R package that takes as input a vertex index of a FreeSurfer brain mesh in MNI305 space (typically fsaverage) and identifies the location in different coordinate systems and with respect to different brain regions. We use this to report the exact locations of clusters or other differences we find. It can also be used to find the vertex in an fsaverage mesh that is closest to a given MNI305 coordinate.

### Coordinate transformation

Supported coordinate systems include:

* MNI305 space RAS coordinates (simply the coordinate of the input vertex index).
* MNI152 space coordinates computed from MNI305 coordinates, with 2 different methods available:
  - using the linear transformation method with the 4x4 FreeSurfer matrix ([section 8 here](https://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems)), or 
  - the more accurate [regfusionr](https://github.com/dfsp-spirit/regfusionr) method.
* Talairach coordinates using [Matthew Brett's transform](http://brainmap.org/training/BrettTransform.html) from MNI152.

#### Validation

If you want to double-check the results of the coordinate transformations, I recommend to use Freeview (comes with FreeSurfer) in combination with the [MNI to Talairach Tool](https://bioimagesuiteweb.github.io/bisweb-manual/tools/mni2tal.html) from Bioimagesuite, which uses the mapping described in Lacadie *et al.*, Neuroimage. 2008 Aug 15; 42(2): 717–725. 

Here is an example for fsaverage vertex 145029:

![Fig1a](./web/fsaverage_vertex_lh_145029.png?raw=true "Vertex 145029 on the left fsaverage surface.")

**Fig. 1a** *Vertex 145029 on the left fsaverage surface (at pink marker). Screenshot from the FreeView application that comes with [FreeSurfer](https://freesurfer.net).* 

![Fig1b](./web/fsaverage_vertex_lh_145029_MNI152_-39_-30_65.png?raw=true "Vertex 145029 on the left fsaverage surface.")

**Fig. 1b** *Location of MNI coordinate 39 -30  65, the result of mapping fsaverage vertex 145029 to MNI152 space. Screenshot from the [MNI - Talairach Tool](https://bioimagesuiteweb.github.io/bisweb-manual/tools/mni2tal.html).*



### Closest brain atlas regions

For a vertex, the package also computes the closest brain regions and the distances to them, based on an atlas. Given a cluster (or more generally, a set of vertices), it can also compute all atlas regions a cluster overlaps with and the percentage as well as absolute number of cluster vertices in the respective atlas regions.

You can use any atlas you like, the three default ones that [come with FreeSurfer](https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation) are:

* The Desikan-Killiany atlas (Desikan *et al.*, 2006. Neuroimage, 31(3):968-80).
* The Destrieux atlas (Destrieux *et al.*, 2010. Neuroimage, 53(1):1-15).
* The DKT40 altas from the [Mindboggle data set](https://mindboggle.info/data.html).

To use a different or custom atlas, just drop the respective annot files for the two hemispheres into the `/label/` directory of your template subject. See [the FreeSurfer documentation on Cortical Parcellations](https://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation) for details on FreeSurfer brain atlases.


### Atlas distances

The following methods are available to compute the distance of a point on a brain surface to a brain atlas region:

* Euclidean distance: point to mean value of region coordinates
* Euclidean distance: point to closest vertex of region
* Geodesic distance to closest vertex of region

