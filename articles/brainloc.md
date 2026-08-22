# brainloc

## About brainloc

The brainloc package serves to describe a location on a brain surface
mesh that is in MNI305 space, like the surfaces of the fsaverage
standard brain templates used by FreeSurfer. Typical use cases of this
package in the field of computational neuroimaging include:

- Finding the mesh vertex closest to a given coordinate (using Euclidean
  distance).
- Finding the coordinates of a vertex on in various standard coordinate
  systems (MNI305, MNI152 and Talairach).
- Finding the brain region the vertex lies in (according to any brain
  atlas available for the respective surface, like the Desikan atlas).
- Finding the distances of all atlas regions to the given vertex (using
  Euclidean or geodesic distance).

Statistical results are often presented as significant clusters of
activation on a brain surface mesh (sets of connected vertices). Typical
use cases for working with clusters are:

- Finding the extreme value of a cluster (min, max, or absolute max) and
  the vertex at which it occurs.
- Finding the peaks within a cluster.
- Finding the brain regions that overlap with the cluster (including the
  degree of overlap).

This vignette describes how to perform these tasks with the brainloc
package.

### Preparation: Getting the data required to run the examples

You will need a directory with the brain templates you want to use,
typically these are the fsaverage subjects that come with the FreeSurfer
neuroimaging software suite. There are 2 ways to get them:

#### Using the template subjects from a local FreeSurfer installation

When FreeSurfer is installed and configured properly on your system,
there is nothing you need to do. In that case, the environment variable
`SUBJECTS_DIR` points at the `subjects` directory under the installation
directory (`FREESURFER_HOME`), and it will be used to find the templates
by the `get_subjects_dir` function. You can also set `SUBJECTS_DIR` to
your custom template directory, of course, and that will be used
instead.

If you cannot or do not want to set environment variables, you can set
the option `brainloc.subjects_dir` in R to any directory that contains
the expected directory structure like this:
`options("brainloc.subjects_dir"="~/my_fs_home/subjects/")`. This would
mean that `~/my_fs_home/subjects/` exists on your system and holds the
directories of the template subjects. We assume that you have the
`fsaverage` subject in there. Setting this option takes precedence over
any environment variables.

#### Have the brainloc package download them for you from the internet

This requires that you have an internet connection and accept the
FreeSurfer software license. See the code in the first example of this
vignette below, which downloads the templates. The templates will be
downloaded only if they cannot be found, and downloaded ones will be
found on the next call to the function, of course.

## Usage Examples Part I: Location information for a vertex or coordinate

### Creating a brainparc object

Let’s say we have a coordinate in MNI305 surface space and want to find
the closest vertex. We first load `brainloc` and create a `brainparc`
object, which internally is a named list holding the surfaces (the left
and right hemisphere meshes) and the surface parcellations (according to
some brain atlas) for these meshes. The parcellations assign to each
vertex of a mesh a single region.

One can construct a brainparc manually from arbitrary loaded meshes and
parcellations using the `brainparc` function, but here we use the
function `brainparc_fs`, which is more convenient if you want to read it
from files in a `SUBJECTS_DIR` of `recon-all` output data. We create a
brainparc for the `fsaverage` template:

``` r

library("brainloc");
# Searches for FreeSurfer SUBJECTS_DIR on machine and downloads if needed.
sjd = get_subjects_dir(allow_download = TRUE, accept_freesurfer_license = TRUE); 
```

    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/label/lh.aparc.a2009s.annot' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/label/lh.aparc.a2009s.annot'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/label/rh.aparc.a2009s.annot' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/label/rh.aparc.a2009s.annot'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/label/lh.aparc.annot' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/label/lh.aparc.annot'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/label/rh.aparc.annot' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/label/rh.aparc.annot'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/label/lh.cortex.label' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/label/lh.cortex.label'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/label/rh.cortex.label' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/label/rh.cortex.label'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/mri/brain.mgz' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/mri/brain.mgz'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/surf/lh.white' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/surf/lh.white'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/surf/rh.white' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/surf/rh.white'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/surf/lh.pial' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/surf/lh.pial'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/surf/rh.pial' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/surf/rh.pial'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/surf/lh.inflated' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/surf/lh.inflated'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/surf/rh.inflated' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/surf/rh.inflated'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/surf/lh.curv' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/surf/lh.curv'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/surf/rh.curv' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/surf/rh.curv'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/ext/FreeSurferColorLUT.txt' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/ext/FreeSurferColorLUT.txt'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage/LICENSE' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage/LICENSE'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage3/label/lh.cortex.label' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage3/label/lh.cortex.label'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage3/label/rh.cortex.label' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage3/label/rh.cortex.label'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage3/surf/lh.white' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage3/surf/lh.white'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage3/surf/rh.white' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage3/surf/rh.white'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/fsaverage3/LICENSE' from 'http://rcmd.org/projects/nitestdata/subjects_dir/fsaverage3/LICENSE'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/subject1/surf/lh.thickness.fwhm0.fsaverage3.mgz' from 'http://rcmd.org/projects/nitestdata/subjects_dir/subject1/surf/lh.thickness.fwhm0.fsaverage3.mgz'
    ## Download file to '/home/runner/.local/share/R/brainloc/subjects_dir/subject1/surf/rh.thickness.fwhm0.fsaverage3.mgz' from 'http://rcmd.org/projects/nitestdata/subjects_dir/subject1/surf/rh.thickness.fwhm0.fsaverage3.mgz'

``` r

bp = brainparc_fs(sjd, "fsaverage", surface = "white", atlas = "aparc");
```

By default, many functions in the `brainloc` package print information
to STDOUT, which is handy for creating command line clients from the
library. We do not want that here, so we disable it globally here:

``` r

options("brainloc.silent"=TRUE);
```

Note: One can also decide for each function that prints something, they
all have a `silent` parameter.

### Finding the vertex closest to a query coordinate

Let’s say we have some query coordinates in surface space and want to
find the closest vertex in the `brainparc` for each of the coordinates.
Here is how we can do that:

``` r

query_coords = matrix(seq(9)+0.1, ncol = 3, byrow = TRUE);
query_coords
```

    ##      [,1] [,2] [,3]
    ## [1,]  1.1  2.1  3.1
    ## [2,]  4.1  5.1  6.1
    ## [3,]  7.1  8.1  9.1

``` r

ccv = coord_closest_vertex(query_coords, get_surface(bp));
ccv
```

    ##   query_x query_y query_z lh_closest_vertex lh_distance rh_closest_vertex
    ## 1     1.1     2.1     3.1            126954   0.8923483            111013
    ## 2     4.1     5.1     6.1             69364   3.7896230             20707
    ## 3     7.1     8.1     9.1             39537   6.7494234             45778
    ##   rh_distance both_closest_vertex both_distance both_hemi
    ## 1   0.4535498              111013     0.4535498        rh
    ## 2   2.9344436               20707     2.9344436        rh
    ## 3   6.2132350               45778     6.2132350        rh

This assumes that the query coordinates are in the surfaces space and
uses Euclidean distance. It works for all subjects and meshes.

### Finding the MNI152 and Talairach coordinates for a point or vertex in MNI305 space

This assumes that you query vertex (or coordinate) is in MNI305 space.
This is true for all FreeSurfer standard template subjects, like
`fsaverage`, `fsaverage6`, etc. Our `brainparc` is for `fsaverage`, so
we can use the coordinates of its vertices to query.

``` r

query_vertices = c(10L, 145029L);
query_hemis = c("lh", "rh"); # vertex 10 on left hemi and 145029L on right hemi.
query_coords = get_surface_coords(bp, query_vertices, query_hemis);
query_coords
```

    ##            [,1]      [,2]      [,3]
    ## [1,] -26.305798  23.81276 -6.026234
    ## [2,]   8.849497 -35.22620 65.073280

``` r

coord_info = coord_MNI305_info(query_coords, surface = get_surface(bp));
coord_info
```

    ## $mni305
    ##            [,1]      [,2]      [,3]
    ## [1,] -26.305798  23.81276 -6.026234
    ## [2,]   8.849497 -35.22620 65.073280
    ## 
    ## $mni152
    ##           [,1]      [,2]      [,3]
    ## [1,] -26.56283  25.01419 -4.704241
    ## [2,]  10.18691 -33.73528 66.281128
    ## 
    ## $talairach
    ##           [,1]      [,2]      [,3]
    ## [1,] -26.29720  24.03707 -5.160009
    ## [2,]  10.08504 -29.63502 62.541625

Note: If you have the `regfusionr` package installed, have a look at the
`method` parameter of the `coord_MNI305_info` function: you can use the
more accurate `regfusionr` method to transform from MNI305 to MNI152
then. Open the docs with
[`?coord_MNI305_info`](https://dfsp-spirit.github.io/brainloc/reference/coord_MNI305_info.md)
to learn more.

### Finding the brain atlas region a vertex belongs to, and the distance to all other brain atlas regions from the vertex

We can easily find the regions a vertex is assigned to in all atlases of
a brain parcellation:

``` r

query_vertices = c(10L, 145029L);
query_hemis = c("lh", "rh"); # vertex 10 on left hemi and 145029L on right hemi.
regions = vertex_regions(bp, query_vertices, query_hemis);

regions
```

    ##   vertex hemi atlas        vertex_region
    ## 1     10   lh aparc lateralorbitofrontal
    ## 2 145029   rh aparc          postcentral

One can also compute the distance of a vertex to all atlas regions.
Several definitions for the distance of a vertex to a region are
available, here we go with the defaults:

``` r

query_vertices = c(10L, 145029L);
query_hemis = c("lh", "rh"); # vertex 10 on left hemi and 145029L on right hemi.
region_dists = vertex_closest_regions(bp, query_vertices, query_hemis, num_regions_to_report=5L);

region_dists
```

    ##    vertex hemi linkage  distance atlas vertex_region_direct vertex_region_nth
    ## 1      10   lh  single euclidean aparc lateralorbitofrontal                 1
    ## 2      10   lh  single euclidean aparc lateralorbitofrontal                 2
    ## 3      10   lh  single euclidean aparc lateralorbitofrontal                 3
    ## 4      10   lh  single euclidean aparc lateralorbitofrontal                 4
    ## 5      10   lh  single euclidean aparc lateralorbitofrontal                 5
    ## 6  145029   rh  single euclidean aparc          postcentral                 1
    ## 7  145029   rh  single euclidean aparc          postcentral                 2
    ## 8  145029   rh  single euclidean aparc          postcentral                 3
    ## 9  145029   rh  single euclidean aparc          postcentral                 4
    ## 10 145029   rh  single euclidean aparc          postcentral                 5
    ##    vertex_region_n_name vertex_region_n_distance vertex_region_n_point_x
    ## 1  lateralorbitofrontal                 0.000000              -26.305798
    ## 2                insula                 5.483304              -29.827978
    ## 3      parstriangularis                10.907174              -34.240818
    ## 4         parsorbitalis                14.363726              -39.584057
    ## 5   medialorbitofrontal                16.396578              -12.834730
    ## 6           postcentral                 0.000000                8.849497
    ## 7            precentral                 1.547366                7.996873
    ## 8           paracentral                 4.240139                9.425453
    ## 9             precuneus                 6.363814               12.051294
    ## 10     superiorparietal                 8.489906               16.539742
    ##    vertex_region_n_point_y vertex_region_n_point_z
    ## 1                 23.81276               -6.026234
    ## 2                 20.04161               -4.171691
    ## 3                 30.83913               -3.450967
    ## 4                 29.04067               -7.661278
    ## 5                 20.36029              -14.712927
    ## 6                -35.22620               65.073280
    ## 7                -33.93639               65.134644
    ## 8                -39.40303               64.624718
    ## 9                -39.72830               68.232040
    ## 10               -37.11887               68.132111

One can change the distance definition with the parameters `distance`
and `linkage` of the `vertex_closest_regions` function. Open the docs
with
[`?vertex_closest_regions`](https://dfsp-spirit.github.io/brainloc/reference/vertex_closest_regions.md)
to learn more.

## Usage Examples Part II: Location information for a cluster

A cluster is a set of vertices on a surface, and the statistical values
that are assigned to the cluster’s vertices. Clusters are used in
vertex-wise analyses to define activated brain regions (where the exact
meaning of *activated* depends on the context). A clusters has a
freeform label or name, but typically only consecutive numbers are used.

The data structure used to store clusters in `brainloc` is called
`clusterinfo`. The following two concepts are relevant to understanding
how a `clusterinfo` instance is typically constructed:

- a `cluster overlay` is an integer vector that assigns to each surface
  vertex a cluster number (zero is assigned to background/non-cluster
  vertices).
- a `statistical map` contains arbitrary values, one for each surface
  vertex (a.k.a. per-vertex data). A *thresholded* map only contains
  values for the cluster vertices, the values for all non-cluster
  vertices are set to an arbitrary value (typically *0.0*).

A `clusterinfo` object can be constructed for a subject in two ways:

- if you have both a `cluster overlay` and the `statistical map`
  (thresholded or unthresholded), use the `clusterinfo` function.
- if you have only a thresholded `statistical map`, use the
  `clusterinfo_from_thresholded_overlay` function instead.

In the following example, we construct a `clusterinfo` instance from a
cluster overlay and the statistical map:

``` r

lh_tmap_file = system.file("extdata", "lh.tmap.mgh", package = "brainloc", mustWork = TRUE);
rh_tmap_file = system.file("extdata", "rh.tmap.mgh", package = "brainloc", mustWork = TRUE);
lh_overlay_file = system.file("extdata", "lh.cluster.overlayID.mgh", package = "brainloc", mustWork = TRUE);
rh_overlay_file = system.file("extdata", "rh.cluster.overlayID.mgh", package = "brainloc", mustWork = TRUE);
clinfo = clusterinfo(lh_overlay_file, rh_overlay_file, lh_tmap_file, rh_tmap_file, template_subject = "fsaverage", subjects_dir = sjd);
is.clusterinfo(clinfo);
```

    ## [1] TRUE

One can extract the clusters from the \`clusterinfo\`\` instance like
this:

``` r

clusters = get_clusters(clinfo);
```

    ## Warning in rgl.init(initValue, onlyNULL): RGL: unable to open X11 display

    ## Warning: 'rgl.init' failed, will use the null device.
    ## See '?rgl.useNULL' for ways to avoid this warning.

This returns a named list, where the keys are cluster names and the
values are integer vectors containing the vertices that belong the the
respective cluster. See
[`?get_clusters`](https://dfsp-spirit.github.io/brainloc/reference/get_clusters.md)
for more options, like returning only the clusters for a specific
hemisphere.

The `clusterinfo` object can now be used to compute various information
on the clusters. We will start by computing information on the location
of the clusters:

``` r

cl_details = cluster_location_details(clinfo);
```

    ## Download file to '/home/runner/.local/share/R/brainloc/talairach/talairach.nii' from 'http://www.talairach.org/talairach.nii'
    ## Download file to '/home/runner/.local/share/R/brainloc/talairach/labels.txt' from 'http://www.talairach.org/labels.txt'

``` r

cl_details;
```

    ##        cluster hemi num_vertices extremum_value extremum_vertex  mni152_r
    ## 1 lh_cluster_2   lh        15889      -3.612464           36395 -56.44855
    ## 2 lh_cluster_1   lh         1949       3.034767          124693 -38.24387
    ##    mni152_a  mni152_s talairach_r talairach_a talairach_s talairach_label_lvl1
    ## 1 -31.87322  40.38953   -55.88407   -29.02201    38.65962                 <NA>
    ## 2  25.68162 -14.75147   -37.86143    24.26172   -13.62192                 <NA>
    ##   talairach_label_lvl2 talairach_label_lvl3 talairach_label_lvl4
    ## 1                 <NA>                 <NA>                 <NA>
    ## 2                 <NA>                 <NA>                 <NA>
    ##   talairach_label_lvl5 talairach_label_full         atlas_region
    ## 1                 <NA>                 <NA>        supramarginal
    ## 2                 <NA>                 <NA> lateralorbitofrontal

Note that the `clusterinfo` object contains a brain parcellation (a
`brainparc` object, see at the top). We can also compute the overlap of
the cluster with its atlas regions:

``` r

overlap = cluster_region_overlap(clinfo);
overlap;
```

    ##                  region num_shared_vertices percent_shared_vertices
    ## 3           postcentral                4880               30.713072
    ## 4            precentral                4224               26.584429
    ## 7         supramarginal                3963               24.941784
    ## 6      superiorparietal                1501                9.446787
    ## 1   caudalmiddlefrontal                 609                3.832840
    ## 2      inferiorparietal                 371                2.334949
    ## 5  rostralmiddlefrontal                 341                2.146139
    ## 21 lateralorbitofrontal                1146               58.799384
    ## 31        parsorbitalis                 629               32.272960
    ## 41     parstriangularis                 115                5.900462
    ## 11               insula                  59                3.027193
    ##    cluster_percent_of_region hemi      cluster atlas
    ## 3                  51.265889   lh lh_cluster_2 aparc
    ## 4                  39.329609   lh lh_cluster_2 aparc
    ## 7                  46.081395   lh lh_cluster_2 aparc
    ## 6                  14.355394   lh lh_cluster_2 aparc
    ## 1                  16.300857   lh lh_cluster_2 aparc
    ## 2                   4.713505   lh lh_cluster_2 aparc
    ## 5                   4.707994   lh lh_cluster_2 aparc
    ## 21                 27.363897   lh lh_cluster_1 aparc
    ## 31                 65.794979   lh lh_cluster_1 aparc
    ## 41                  5.620723   lh lh_cluster_1 aparc
    ## 11                  1.128323   lh lh_cluster_1 aparc

Sometimes one wants to find all peaks of a cluster. A peak exists at
vertices which have the most extreme value in their neighborhood. For a
brain mesh, the neighborhood of a vertex is defined as all vertices
which are reachable by one hop along the mesh edges.

``` r

peaks = cluster_peaks(clinfo);
peaks
```

    ##          cluster hemi peak_vertex peak_value
    ## 1   lh_cluster_2   lh          17  -2.744687
    ## 2   lh_cluster_2   lh          48  -1.858908
    ## 3   lh_cluster_2   lh         167  -3.017634
    ## 4   lh_cluster_2   lh         188  -2.706193
    ## 5   lh_cluster_2   lh         368  -3.061095
    ## 6   lh_cluster_2   lh         844  -2.856204
    ## 7   lh_cluster_2   lh        1059  -2.536362
    ## 8   lh_cluster_2   lh        1443  -2.452701
    ## 9   lh_cluster_2   lh        1570  -2.675941
    ## 10  lh_cluster_2   lh        1818  -3.056978
    ## 11  lh_cluster_2   lh        2096  -2.821692
    ## 12  lh_cluster_2   lh        2633  -2.942444
    ## 13  lh_cluster_2   lh        2939  -3.173836
    ## 14  lh_cluster_2   lh        2946  -2.900167
    ## 15  lh_cluster_2   lh        3368  -1.801927
    ## 16  lh_cluster_2   lh        3382  -2.089746
    ## 17  lh_cluster_2   lh        3445  -2.607287
    ## 18  lh_cluster_2   lh        4006  -2.744771
    ## 19  lh_cluster_2   lh        4482  -2.015952
    ## 20  lh_cluster_2   lh        4507  -2.937274
    ## 21  lh_cluster_2   lh        4604  -2.673632
    ## 22  lh_cluster_2   lh        4609  -2.855095
    ## 23  lh_cluster_2   lh        5774  -2.092039
    ## 24  lh_cluster_2   lh        6137  -2.874391
    ## 25  lh_cluster_2   lh        6138  -2.841173
    ## 26  lh_cluster_2   lh        6140  -2.883987
    ## 27  lh_cluster_2   lh        6145  -2.785391
    ## 28  lh_cluster_2   lh        6446  -2.985426
    ## 29  lh_cluster_2   lh        7220  -2.988027
    ## 30  lh_cluster_2   lh        7227  -2.464265
    ## 31  lh_cluster_2   lh        7492  -2.424789
    ## 32  lh_cluster_2   lh        7504  -2.538386
    ## 33  lh_cluster_2   lh        7695  -3.172078
    ## 34  lh_cluster_2   lh        8081  -3.030712
    ## 35  lh_cluster_2   lh        8115  -2.020023
    ## 36  lh_cluster_2   lh        8125  -2.621741
    ## 37  lh_cluster_2   lh        8151  -1.715064
    ## 38  lh_cluster_2   lh        8733  -3.057757
    ## 39  lh_cluster_2   lh        8827  -2.437251
    ## 40  lh_cluster_2   lh        8834  -1.913942
    ## 41  lh_cluster_2   lh        8836  -2.157695
    ## 42  lh_cluster_2   lh        8854  -2.799692
    ## 43  lh_cluster_2   lh        8858  -2.493563
    ## 44  lh_cluster_2   lh       10358  -2.780647
    ## 45  lh_cluster_2   lh       10454  -2.474493
    ## 46  lh_cluster_2   lh       10619  -3.218110
    ## 47  lh_cluster_2   lh       10647  -1.899229
    ## 48  lh_cluster_2   lh       11048  -1.802100
    ## 49  lh_cluster_2   lh       11071  -2.285401
    ## 50  lh_cluster_2   lh       11123  -2.745003
    ## 51  lh_cluster_2   lh       11671  -2.953016
    ## 52  lh_cluster_2   lh       11731  -1.866850
    ## 53  lh_cluster_2   lh       11890  -2.089939
    ## 54  lh_cluster_2   lh       12287  -2.859982
    ## 55  lh_cluster_2   lh       13449  -1.874903
    ## 56  lh_cluster_2   lh       13528  -2.266294
    ## 57  lh_cluster_2   lh       13616  -2.872344
    ## 58  lh_cluster_2   lh       13656  -2.923451
    ## 59  lh_cluster_2   lh       14073  -2.930758
    ## 60  lh_cluster_2   lh       14124  -2.218930
    ## 61  lh_cluster_2   lh       14580  -3.094419
    ## 62  lh_cluster_2   lh       14588  -3.006916
    ## 63  lh_cluster_2   lh       15891  -2.470590
    ## 64  lh_cluster_2   lh       15893  -2.477159
    ## 65  lh_cluster_2   lh       16094  -2.643576
    ## 66  lh_cluster_2   lh       16330  -2.195119
    ## 67  lh_cluster_2   lh       16381  -3.135160
    ## 68  lh_cluster_2   lh       16386  -2.641017
    ## 69  lh_cluster_2   lh       16468  -2.387602
    ## 70  lh_cluster_2   lh       17688  -1.885616
    ## 71  lh_cluster_2   lh       17729  -3.053454
    ## 72  lh_cluster_2   lh       17731  -3.057509
    ## 73  lh_cluster_2   lh       17735  -3.044742
    ## 74  lh_cluster_2   lh       17739  -3.058966
    ## 75  lh_cluster_2   lh       17781  -1.703062
    ## 76  lh_cluster_2   lh       17798  -1.951260
    ## 77  lh_cluster_2   lh       17819  -1.723828
    ## 78  lh_cluster_2   lh       17877  -1.749732
    ## 79  lh_cluster_2   lh       17878  -1.769274
    ## 80  lh_cluster_2   lh       19085  -1.838858
    ## 81  lh_cluster_2   lh       19139  -2.318835
    ## 82  lh_cluster_2   lh       19150  -1.908982
    ## 83  lh_cluster_2   lh       19242  -3.534949
    ## 84  lh_cluster_2   lh       19250  -2.790674
    ## 85  lh_cluster_2   lh       19266  -2.932451
    ## 86  lh_cluster_2   lh       19267  -2.913095
    ## 87  lh_cluster_2   lh       19278  -2.867075
    ## 88  lh_cluster_2   lh       19317  -1.953748
    ## 89  lh_cluster_2   lh       19517  -2.446122
    ## 90  lh_cluster_2   lh       19583  -2.817948
    ## 91  lh_cluster_2   lh       19599  -2.159396
    ## 92  lh_cluster_2   lh       23052  -1.897364
    ## 93  lh_cluster_2   lh       23448  -1.892832
    ## 94  lh_cluster_2   lh       23552  -2.558031
    ## 95  lh_cluster_2   lh       23689  -1.806096
    ## 96  lh_cluster_2   lh       23719  -3.149762
    ## 97  lh_cluster_2   lh       23720  -3.107920
    ## 98  lh_cluster_2   lh       23725  -3.252490
    ## 99  lh_cluster_2   lh       23727  -3.067868
    ## 100 lh_cluster_2   lh       24488  -1.851189
    ## 101 lh_cluster_2   lh       24490  -1.823261
    ## 102 lh_cluster_2   lh       24504  -2.965343
    ## 103 lh_cluster_2   lh       24518  -2.093884
    ## 104 lh_cluster_2   lh       24777  -1.797602
    ## 105 lh_cluster_2   lh       24850  -2.212785
    ## 106 lh_cluster_2   lh       24853  -2.454756
    ## 107 lh_cluster_2   lh       24858  -2.361548
    ## 108 lh_cluster_2   lh       24866  -2.687023
    ## 109 lh_cluster_2   lh       25569  -2.793199
    ## 110 lh_cluster_2   lh       25609  -2.821468
    ## 111 lh_cluster_2   lh       25647  -1.786969
    ## 112 lh_cluster_2   lh       25714  -1.688275
    ## 113 lh_cluster_2   lh       26519  -2.699700
    ## 114 lh_cluster_2   lh       26542  -2.936082
    ## 115 lh_cluster_2   lh       26569  -1.858159
    ## 116 lh_cluster_2   lh       26711  -2.694932
    ## 117 lh_cluster_2   lh       28822  -2.944403
    ## 118 lh_cluster_2   lh       28840  -1.971871
    ## 119 lh_cluster_2   lh       29098  -2.386630
    ## 120 lh_cluster_2   lh       29299  -2.876646
    ## 121 lh_cluster_2   lh       29680  -2.922689
    ## 122 lh_cluster_2   lh       29698  -2.660622
    ## 123 lh_cluster_2   lh       30404  -2.991887
    ## 124 lh_cluster_2   lh       30406  -3.056174
    ## 125 lh_cluster_2   lh       30537  -2.745333
    ## 126 lh_cluster_2   lh       32177  -1.789202
    ## 127 lh_cluster_2   lh       32178  -1.733267
    ## 128 lh_cluster_2   lh       32197  -2.825445
    ## 129 lh_cluster_2   lh       32198  -2.897207
    ## 130 lh_cluster_2   lh       32210  -3.026260
    ## 131 lh_cluster_2   lh       32480  -2.676601
    ## 132 lh_cluster_2   lh       32490  -2.677931
    ## 133 lh_cluster_2   lh       32491  -2.654871
    ## 134 lh_cluster_2   lh       32765  -1.849113
    ## 135 lh_cluster_2   lh       32785  -1.838812
    ## 136 lh_cluster_2   lh       32823  -2.453201
    ## 137 lh_cluster_2   lh       32830  -2.649307
    ## 138 lh_cluster_2   lh       32966  -2.562637
    ## 139 lh_cluster_2   lh       34471  -1.865558
    ## 140 lh_cluster_2   lh       34574  -3.084813
    ## 141 lh_cluster_2   lh       34583  -2.798079
    ## 142 lh_cluster_2   lh       34634  -1.920266
    ## 143 lh_cluster_2   lh       34648  -1.902689
    ## 144 lh_cluster_2   lh       34651  -1.949643
    ## 145 lh_cluster_2   lh       34654  -1.997008
    ## 146 lh_cluster_2   lh       34657  -1.956390
    ## 147 lh_cluster_2   lh       34663  -1.997914
    ## 148 lh_cluster_2   lh       34666  -2.019443
    ## 149 lh_cluster_2   lh       34723  -2.339954
    ## 150 lh_cluster_2   lh       34724  -2.279819
    ## 151 lh_cluster_2   lh       34725  -2.269535
    ## 152 lh_cluster_2   lh       34756  -1.704730
    ## 153 lh_cluster_2   lh       36203  -1.887507
    ## 154 lh_cluster_2   lh       36278  -2.215168
    ## 155 lh_cluster_2   lh       36346  -1.881644
    ## 156 lh_cluster_2   lh       36350  -1.963081
    ## 157 lh_cluster_2   lh       36395  -3.612464
    ## 158 lh_cluster_2   lh       36404  -2.869724
    ## 159 lh_cluster_2   lh       36418  -2.742117
    ## 160 lh_cluster_2   lh       36425  -2.865510
    ## 161 lh_cluster_2   lh       36433  -2.899261
    ## 162 lh_cluster_2   lh       36441  -3.034779
    ## 163 lh_cluster_2   lh       36443  -3.019564
    ## 164 lh_cluster_2   lh       36794  -2.935117
    ## 165 lh_cluster_2   lh       36814  -2.776366
    ## 166 lh_cluster_2   lh       36820  -2.459363
    ## 167 lh_cluster_2   lh       41033  -2.907137
    ## 168 lh_cluster_2   lh       41163  -2.492348
    ## 169 lh_cluster_2   lh       41179  -1.820250
    ## 170 lh_cluster_2   lh       41196  -1.685994
    ## 171 lh_cluster_2   lh       41323  -1.865373
    ## 172 lh_cluster_2   lh       41338  -2.323637
    ## 173 lh_cluster_2   lh       41341  -3.178528
    ## 174 lh_cluster_2   lh       41788  -2.340587
    ## 175 lh_cluster_2   lh       41793  -1.786912
    ## 176 lh_cluster_2   lh       41936  -2.769045
    ## 177 lh_cluster_2   lh       41974  -2.349726
    ## 178 lh_cluster_2   lh       42451  -1.843482
    ## 179 lh_cluster_2   lh       42925  -1.927888
    ## 180 lh_cluster_2   lh       44378  -2.907199
    ## 181 lh_cluster_2   lh       44768  -2.641192
    ## 182 lh_cluster_2   lh       44791  -2.997261
    ## 183 lh_cluster_2   lh       44839  -1.948231
    ## 184 lh_cluster_2   lh       44843  -2.225387
    ## 185 lh_cluster_2   lh       44964  -2.287409
    ## 186 lh_cluster_2   lh       45288  -2.194043
    ## 187 lh_cluster_2   lh       45304  -3.026557
    ## 188 lh_cluster_2   lh       45413  -2.396898
    ## 189 lh_cluster_2   lh       45824  -2.435784
    ## 190 lh_cluster_2   lh       46617  -2.202813
    ## 191 lh_cluster_2   lh       46815  -2.667459
    ## 192 lh_cluster_2   lh       46821  -2.445037
    ## 193 lh_cluster_2   lh       47173  -2.010066
    ## 194 lh_cluster_2   lh       47182  -2.365768
    ## 195 lh_cluster_2   lh       47185  -2.387568
    ## 196 lh_cluster_2   lh       47209  -2.548884
    ## 197 lh_cluster_2   lh       48372  -2.871494
    ## 198 lh_cluster_2   lh       48385  -1.697148
    ## 199 lh_cluster_2   lh       48403  -1.791297
    ## 200 lh_cluster_2   lh       48406  -1.788480
    ## 201 lh_cluster_2   lh       48457  -3.074100
    ## 202 lh_cluster_2   lh       48461  -3.084569
    ## 203 lh_cluster_2   lh       48477  -2.853503
    ## 204 lh_cluster_2   lh       48486  -2.693330
    ## 205 lh_cluster_2   lh       48516  -1.975662
    ## 206 lh_cluster_2   lh       48520  -1.947287
    ## 207 lh_cluster_2   lh       48557  -2.614529
    ## 208 lh_cluster_2   lh       48619  -1.766548
    ## 209 lh_cluster_2   lh       49945  -2.227815
    ## 210 lh_cluster_2   lh       49964  -3.046242
    ## 211 lh_cluster_2   lh       49970  -2.802089
    ## 212 lh_cluster_2   lh       49984  -2.960794
    ## 213 lh_cluster_2   lh       49992  -2.972959
    ## 214 lh_cluster_2   lh       50211  -1.840944
    ## 215 lh_cluster_2   lh       50255  -1.891566
    ## 216 lh_cluster_2   lh       50284  -2.756488
    ## 217 lh_cluster_2   lh       50287  -2.889412
    ## 218 lh_cluster_2   lh       50291  -2.891750
    ## 219 lh_cluster_2   lh       50295  -2.901865
    ## 220 lh_cluster_2   lh       53977  -2.900434
    ## 221 lh_cluster_2   lh       54387  -2.997250
    ## 222 lh_cluster_2   lh       54391  -3.097589
    ## 223 lh_cluster_2   lh       54406  -2.460313
    ## 224 lh_cluster_2   lh       54849  -2.191310
    ## 225 lh_cluster_2   lh       54852  -2.165772
    ## 226 lh_cluster_2   lh       54859  -2.085207
    ## 227 lh_cluster_2   lh       54892  -3.203849
    ## 228 lh_cluster_2   lh       55001  -2.402891
    ## 229 lh_cluster_2   lh       55421  -3.095200
    ## 230 lh_cluster_2   lh       56210  -2.347468
    ## 231 lh_cluster_2   lh       56212  -2.423478
    ## 232 lh_cluster_2   lh       56217  -2.070493
    ## 233 lh_cluster_2   lh       56247  -2.266964
    ## 234 lh_cluster_2   lh       56254  -1.685358
    ## 235 lh_cluster_2   lh       56417  -2.353697
    ## 236 lh_cluster_2   lh       56688  -2.654236
    ## 237 lh_cluster_2   lh       56795  -2.407562
    ## 238 lh_cluster_2   lh       56808  -2.691053
    ## 239 lh_cluster_2   lh       58032  -3.183473
    ## 240 lh_cluster_2   lh       58041  -2.813020
    ## 241 lh_cluster_2   lh       58088  -2.749175
    ## 242 lh_cluster_2   lh       58090  -2.756687
    ## 243 lh_cluster_2   lh       58100  -1.808172
    ## 244 lh_cluster_2   lh       58114  -1.724980
    ## 245 lh_cluster_2   lh       58127  -1.980369
    ## 246 lh_cluster_2   lh       58163  -3.174872
    ## 247 lh_cluster_2   lh       58177  -1.917051
    ## 248 lh_cluster_2   lh       59548  -2.886985
    ## 249 lh_cluster_2   lh       59552  -2.833109
    ## 250 lh_cluster_2   lh       59579  -2.332160
    ## 251 lh_cluster_2   lh       59591  -2.930279
    ## 252 lh_cluster_2   lh       59602  -2.837378
    ## 253 lh_cluster_2   lh       59638  -1.821030
    ## 254 lh_cluster_2   lh       59856  -2.141293
    ## 255 lh_cluster_2   lh       59860  -2.176594
    ## 256 lh_cluster_2   lh       59881  -2.743769
    ## 257 lh_cluster_2   lh       59899  -2.887746
    ## 258 lh_cluster_2   lh       59901  -2.810944
    ## 259 lh_cluster_2   lh       59905  -2.797544
    ## 260 lh_cluster_2   lh       59919  -2.180943
    ## 261 lh_cluster_2   lh       61134  -1.871940
    ## 262 lh_cluster_2   lh       63383  -1.874838
    ## 263 lh_cluster_2   lh       63408  -2.280779
    ## 264 lh_cluster_2   lh       63416  -1.811714
    ## 265 lh_cluster_2   lh       63612  -2.459620
    ## 266 lh_cluster_2   lh       63615  -2.382020
    ## 267 lh_cluster_2   lh       63617  -2.290917
    ## 268 lh_cluster_2   lh       63623  -2.362006
    ## 269 lh_cluster_2   lh       64013  -2.796220
    ## 270 lh_cluster_2   lh       64420  -2.829930
    ## 271 lh_cluster_2   lh       65174  -2.698698
    ## 272 lh_cluster_2   lh       65178  -2.838101
    ## 273 lh_cluster_2   lh       65235  -2.851224
    ## 274 lh_cluster_2   lh       65238  -2.828228
    ## 275 lh_cluster_2   lh       65247  -2.846072
    ## 276 lh_cluster_2   lh       65264  -2.748195
    ## 277 lh_cluster_2   lh       65281  -2.638127
    ## 278 lh_cluster_2   lh       65394  -1.939644
    ## 279 lh_cluster_2   lh       65425  -1.819796
    ## 280 lh_cluster_2   lh       66589  -1.731564
    ## 281 lh_cluster_2   lh       66742  -3.090036
    ## 282 lh_cluster_2   lh       66771  -3.178936
    ## 283 lh_cluster_2   lh       66794  -3.026984
    ## 284 lh_cluster_2   lh       66846  -1.904422
    ## 285 lh_cluster_2   lh       66853  -1.732712
    ## 286 lh_cluster_2   lh       67082  -2.787355
    ## 287 lh_cluster_2   lh       67091  -2.757247
    ## 288 lh_cluster_2   lh       67116  -2.352031
    ## 289 lh_cluster_2   lh       70565  -2.437227
    ## 290 lh_cluster_2   lh       70651  -3.009846
    ## 291 lh_cluster_2   lh       70670  -2.959379
    ## 292 lh_cluster_2   lh       70676  -2.989444
    ## 293 lh_cluster_2   lh       70688  -2.555965
    ## 294 lh_cluster_2   lh       70694  -2.432862
    ## 295 lh_cluster_2   lh       70706  -2.435696
    ## 296 lh_cluster_2   lh       70721  -2.294841
    ## 297 lh_cluster_2   lh       70725  -2.313073
    ## 298 lh_cluster_2   lh       70730  -2.136467
    ## 299 lh_cluster_2   lh       70778  -2.259028
    ## 300 lh_cluster_2   lh       70791  -1.949017
    ## 301 lh_cluster_2   lh       70804  -2.040392
    ## 302 lh_cluster_2   lh       70808  -1.991043
    ## 303 lh_cluster_2   lh       70811  -1.964732
    ## 304 lh_cluster_2   lh       70817  -2.005408
    ## 305 lh_cluster_2   lh       70828  -1.927297
    ## 306 lh_cluster_2   lh       71314  -2.637314
    ## 307 lh_cluster_2   lh       71335  -2.484150
    ## 308 lh_cluster_2   lh       72025  -2.146111
    ## 309 lh_cluster_2   lh       72125  -2.508476
    ## 310 lh_cluster_2   lh       72127  -2.567078
    ## 311 lh_cluster_2   lh       72157  -2.887072
    ## 312 lh_cluster_2   lh       72168  -3.080318
    ## 313 lh_cluster_2   lh       72173  -3.111388
    ## 314 lh_cluster_2   lh       72177  -3.158182
    ## 315 lh_cluster_2   lh       72179  -3.140754
    ## 316 lh_cluster_2   lh       72421  -2.413467
    ## 317 lh_cluster_2   lh       72436  -2.388741
    ## 318 lh_cluster_2   lh       72443  -2.383343
    ## 319 lh_cluster_2   lh       72495  -2.526154
    ## 320 lh_cluster_2   lh       72512  -2.444963
    ## 321 lh_cluster_2   lh       72514  -2.411968
    ## 322 lh_cluster_2   lh       72517  -2.363472
    ## 323 lh_cluster_2   lh       73785  -1.696124
    ## 324 lh_cluster_2   lh       76024  -1.769008
    ## 325 lh_cluster_2   lh       76037  -2.547835
    ## 326 lh_cluster_2   lh       76082  -1.780614
    ## 327 lh_cluster_2   lh       76086  -1.769181
    ## 328 lh_cluster_2   lh       76094  -1.894956
    ## 329 lh_cluster_2   lh       76099  -1.839349
    ## 330 lh_cluster_2   lh       76103  -1.830023
    ## 331 lh_cluster_2   lh       76169  -3.044477
    ## 332 lh_cluster_2   lh       76174  -3.167114
    ## 333 lh_cluster_2   lh       76206  -2.634015
    ## 334 lh_cluster_2   lh       76219  -2.967653
    ## 335 lh_cluster_2   lh       76225  -3.048519
    ## 336 lh_cluster_2   lh       76236  -3.072217
    ## 337 lh_cluster_2   lh       76239  -3.021285
    ## 338 lh_cluster_2   lh       76243  -3.037441
    ## 339 lh_cluster_2   lh       76247  -3.024743
    ## 340 lh_cluster_2   lh       76274  -3.149360
    ## 341 lh_cluster_2   lh       76276  -3.131381
    ## 342 lh_cluster_2   lh       76309  -2.878780
    ## 343 lh_cluster_2   lh       76330  -2.828759
    ## 344 lh_cluster_2   lh       76333  -2.742433
    ## 345 lh_cluster_2   lh       76338  -2.779932
    ## 346 lh_cluster_2   lh       76434  -2.011052
    ## 347 lh_cluster_2   lh       76463  -1.871284
    ## 348 lh_cluster_2   lh       76477  -2.036463
    ## 349 lh_cluster_2   lh       76483  -2.088137
    ## 350 lh_cluster_2   lh       76553  -2.786468
    ## 351 lh_cluster_2   lh       76648  -2.286645
    ## 352 lh_cluster_2   lh       76673  -1.730712
    ## 353 lh_cluster_2   lh       76676  -1.780879
    ## 354 lh_cluster_2   lh       76679  -1.696896
    ## 355 lh_cluster_2   lh       78213  -2.042754
    ## 356 lh_cluster_2   lh       80699  -2.363550
    ## 357 lh_cluster_2   lh       80704  -2.246228
    ## 358 lh_cluster_2   lh       80821  -2.714560
    ## 359 lh_cluster_2   lh       80832  -2.923371
    ## 360 lh_cluster_2   lh       80880  -3.010818
    ## 361 lh_cluster_2   lh       80884  -3.066441
    ## 362 lh_cluster_2   lh       80890  -3.093800
    ## 363 lh_cluster_2   lh       80961  -2.954245
    ## 364 lh_cluster_2   lh       80983  -1.965804
    ## 365 lh_cluster_2   lh       81028  -2.685546
    ## 366 lh_cluster_2   lh       81040  -3.157718
    ## 367 lh_cluster_2   lh       81046  -3.018775
    ## 368 lh_cluster_2   lh       81073  -2.769447
    ## 369 lh_cluster_2   lh       81499  -1.925347
    ## 370 lh_cluster_2   lh       81600  -2.351384
    ## 371 lh_cluster_2   lh       81606  -2.422964
    ## 372 lh_cluster_2   lh       81661  -2.246440
    ## 373 lh_cluster_2   lh       81675  -2.187205
    ## 374 lh_cluster_2   lh       81678  -2.187558
    ## 375 lh_cluster_2   lh       81749  -2.927985
    ## 376 lh_cluster_2   lh       81753  -2.928776
    ## 377 lh_cluster_2   lh       81780  -2.846657
    ## 378 lh_cluster_2   lh       81791  -2.690152
    ## 379 lh_cluster_2   lh       92287  -2.932397
    ## 380 lh_cluster_2   lh       92373  -2.796881
    ## 381 lh_cluster_2   lh       92540  -3.105558
    ## 382 lh_cluster_2   lh       92549  -2.483843
    ## 383 lh_cluster_2   lh       92551  -1.762835
    ## 384 lh_cluster_2   lh       92564  -1.674167
    ## 385 lh_cluster_2   lh       92572  -1.939513
    ## 386 lh_cluster_2   lh       92586  -1.694452
    ## 387 lh_cluster_2   lh       92809  -1.824051
    ## 388 lh_cluster_2   lh       92814  -2.193547
    ## 389 lh_cluster_2   lh       92819  -2.118423
    ## 390 lh_cluster_2   lh       92834  -2.305150
    ## 391 lh_cluster_2   lh       92837  -2.297570
    ## 392 lh_cluster_2   lh       92844  -3.235487
    ## 393 lh_cluster_2   lh       92849  -2.960643
    ## 394 lh_cluster_2   lh       92910  -2.380729
    ## 395 lh_cluster_2   lh       93615  -1.761233
    ## 396 lh_cluster_2   lh       93620  -2.513124
    ## 397 lh_cluster_2   lh       93633  -2.360366
    ## 398 lh_cluster_2   lh       93650  -2.320643
    ## 399 lh_cluster_2   lh       93652  -2.298004
    ## 400 lh_cluster_2   lh       93987  -2.746527
    ## 401 lh_cluster_2   lh       94239  -2.747372
    ## 402 lh_cluster_2   lh       94684  -2.761430
    ## 403 lh_cluster_2   lh       94687  -2.808855
    ## 404 lh_cluster_2   lh       94724  -3.180370
    ## 405 lh_cluster_2   lh       94740  -3.035608
    ## 406 lh_cluster_2   lh       94742  -3.021680
    ## 407 lh_cluster_2   lh       94743  -2.810092
    ## 408 lh_cluster_2   lh       94762  -2.739837
    ## 409 lh_cluster_2   lh       94764  -1.879127
    ## 410 lh_cluster_2   lh       94767  -1.874810
    ## 411 lh_cluster_2   lh       94812  -1.981002
    ## 412 lh_cluster_2   lh       94825  -2.084657
    ## 413 lh_cluster_2   lh       94837  -1.680597
    ## 414 lh_cluster_2   lh       95582  -2.224030
    ## 415 lh_cluster_2   lh       95584  -2.178457
    ## 416 lh_cluster_2   lh       95635  -2.873177
    ## 417 lh_cluster_2   lh       95640  -2.910546
    ## 418 lh_cluster_2   lh       95666  -2.860042
    ## 419 lh_cluster_2   lh       95815  -2.096956
    ## 420 lh_cluster_2   lh       95817  -2.097326
    ## 421 lh_cluster_2   lh       95831  -2.740822
    ## 422 lh_cluster_2   lh       95834  -2.745071
    ## 423 lh_cluster_2   lh       96586  -1.840730
    ## 424 lh_cluster_2   lh       97970  -1.943194
    ## 425 lh_cluster_2   lh       98071  -2.440693
    ## 426 lh_cluster_2   lh       98290  -2.719596
    ## 427 lh_cluster_2   lh       98312  -2.847118
    ## 428 lh_cluster_2   lh       98559  -2.901016
    ## 429 lh_cluster_2   lh       99015  -2.607474
    ## 430 lh_cluster_2   lh       99147  -2.208974
    ## 431 lh_cluster_2   lh       99901  -2.268416
    ## 432 lh_cluster_2   lh       99951  -3.079802
    ## 433 lh_cluster_2   lh       99980  -3.052064
    ## 434 lh_cluster_2   lh      100007  -1.909146
    ## 435 lh_cluster_2   lh      100015  -1.764204
    ## 436 lh_cluster_2   lh      100154  -2.679595
    ## 437 lh_cluster_2   lh      100159  -2.717332
    ## 438 lh_cluster_2   lh      100171  -2.452881
    ## 439 lh_cluster_2   lh      102246  -2.415736
    ## 440 lh_cluster_2   lh      102255  -1.743171
    ## 441 lh_cluster_2   lh      102304  -2.963894
    ## 442 lh_cluster_2   lh      102321  -2.463087
    ## 443 lh_cluster_2   lh      102374  -2.247676
    ## 444 lh_cluster_2   lh      102388  -1.991617
    ## 445 lh_cluster_2   lh      102390  -2.016874
    ## 446 lh_cluster_2   lh      102699  -2.666774
    ## 447 lh_cluster_2   lh      102709  -2.584761
    ## 448 lh_cluster_2   lh      103182  -2.588064
    ## 449 lh_cluster_2   lh      103205  -3.009114
    ## 450 lh_cluster_2   lh      103341  -1.962196
    ## 451 lh_cluster_2   lh      103403  -2.537879
    ## 452 lh_cluster_2   lh      103413  -2.440415
    ## 453 lh_cluster_2   lh      105493  -2.507114
    ## 454 lh_cluster_2   lh      105497  -2.846932
    ## 455 lh_cluster_2   lh      105501  -2.847695
    ## 456 lh_cluster_2   lh      105522  -1.707243
    ## 457 lh_cluster_2   lh      105527  -2.472236
    ## 458 lh_cluster_2   lh      105621  -2.739356
    ## 459 lh_cluster_2   lh      105650  -3.009020
    ## 460 lh_cluster_2   lh      105652  -3.056979
    ## 461 lh_cluster_2   lh      105674  -3.099640
    ## 462 lh_cluster_2   lh      105685  -2.867559
    ## 463 lh_cluster_2   lh      105687  -2.818274
    ## 464 lh_cluster_2   lh      105691  -2.852376
    ## 465 lh_cluster_2   lh      105765  -1.978059
    ## 466 lh_cluster_2   lh      105792  -2.044348
    ## 467 lh_cluster_2   lh      105797  -2.059854
    ## 468 lh_cluster_2   lh      105818  -2.087981
    ## 469 lh_cluster_2   lh      105913  -1.687425
    ## 470 lh_cluster_2   lh      106858  -2.672052
    ## 471 lh_cluster_2   lh      108321  -2.403574
    ## 472 lh_cluster_2   lh      108391  -2.390295
    ## 473 lh_cluster_2   lh      108406  -2.930959
    ## 474 lh_cluster_2   lh      108411  -2.922242
    ## 475 lh_cluster_2   lh      108418  -2.931588
    ## 476 lh_cluster_2   lh      108502  -1.746760
    ## 477 lh_cluster_2   lh      108884  -1.775076
    ## 478 lh_cluster_2   lh      108941  -2.783559
    ## 479 lh_cluster_2   lh      108987  -2.620885
    ## 480 lh_cluster_2   lh      108995  -2.582790
    ## 481 lh_cluster_2   lh      115232  -1.768347
    ## 482 lh_cluster_2   lh      115237  -1.911130
    ## 483 lh_cluster_2   lh      115241  -1.895569
    ## 484 lh_cluster_2   lh      115322  -2.208857
    ## 485 lh_cluster_2   lh      115695  -2.877914
    ## 486 lh_cluster_2   lh      116097  -2.672016
    ## 487 lh_cluster_2   lh      116149  -1.943913
    ## 488 lh_cluster_2   lh      116158  -2.130025
    ## 489 lh_cluster_2   lh      116722  -2.200630
    ## 490 lh_cluster_2   lh      116732  -2.121809
    ## 491 lh_cluster_2   lh      116745  -2.162316
    ## 492 lh_cluster_2   lh      116770  -2.323831
    ## 493 lh_cluster_2   lh      116774  -2.259454
    ## 494 lh_cluster_2   lh      116796  -3.203082
    ## 495 lh_cluster_2   lh      116802  -3.045920
    ## 496 lh_cluster_2   lh      116951  -2.428258
    ## 497 lh_cluster_2   lh      118581  -1.805043
    ## 498 lh_cluster_2   lh      118595  -2.629833
    ## 499 lh_cluster_2   lh      118599  -2.658152
    ## 500 lh_cluster_2   lh      118613  -2.962621
    ## 501 lh_cluster_2   lh      118627  -2.278138
    ## 502 lh_cluster_2   lh      118677  -2.203577
    ## 503 lh_cluster_2   lh      118681  -2.222703
    ## 504 lh_cluster_2   lh      118684  -2.245303
    ## 505 lh_cluster_2   lh      118692  -1.702453
    ## 506 lh_cluster_2   lh      118915  -2.646786
    ## 507 lh_cluster_2   lh      118923  -2.319460
    ## 508 lh_cluster_2   lh      118925  -2.369248
    ## 509 lh_cluster_2   lh      119225  -1.744432
    ## 510 lh_cluster_2   lh      119424  -2.345010
    ## 511 lh_cluster_2   lh      119439  -2.376050
    ## 512 lh_cluster_2   lh      119471  -2.692400
    ## 513 lh_cluster_2   lh      121179  -3.151238
    ## 514 lh_cluster_2   lh      121183  -3.215775
    ## 515 lh_cluster_2   lh      121235  -2.945282
    ## 516 lh_cluster_2   lh      121243  -2.884675
    ## 517 lh_cluster_2   lh      121250  -2.856719
    ## 518 lh_cluster_2   lh      121276  -1.781689
    ## 519 lh_cluster_2   lh      121310  -2.107559
    ## 520 lh_cluster_2   lh      121434  -1.746353
    ## 521 lh_cluster_2   lh      121441  -1.739218
    ## 522 lh_cluster_2   lh      121443  -1.713047
    ## 523 lh_cluster_2   lh      121449  -1.688209
    ## 524 lh_cluster_2   lh      123102  -1.981119
    ## 525 lh_cluster_2   lh      123178  -2.245620
    ## 526 lh_cluster_2   lh      123185  -2.294299
    ## 527 lh_cluster_2   lh      123316  -2.870781
    ## 528 lh_cluster_2   lh      123326  -3.123764
    ## 529 lh_cluster_2   lh      123331  -2.795111
    ## 530 lh_cluster_2   lh      123358  -2.925398
    ## 531 lh_cluster_2   lh      123362  -2.925673
    ## 532 lh_cluster_2   lh      123379  -2.851797
    ## 533 lh_cluster_2   lh      123426  -1.877161
    ## 534 lh_cluster_2   lh      123468  -2.076707
    ## 535 lh_cluster_2   lh      123743  -2.443603
    ## 536 lh_cluster_2   lh      123767  -2.732703
    ## 537 lh_cluster_2   lh      123785  -2.897136
    ## 538 lh_cluster_2   lh      123791  -2.872825
    ## 539 lh_cluster_2   lh      123808  -2.372777
    ## 540 lh_cluster_2   lh      125527  -1.814749
    ## 541 lh_cluster_2   lh      128669  -1.834998
    ## 542 lh_cluster_2   lh      128702  -2.308194
    ## 543 lh_cluster_2   lh      128740  -1.884908
    ## 544 lh_cluster_2   lh      128754  -2.332346
    ## 545 lh_cluster_2   lh      128974  -2.406472
    ## 546 lh_cluster_2   lh      129418  -2.699180
    ## 547 lh_cluster_2   lh      129425  -2.821688
    ## 548 lh_cluster_2   lh      129454  -2.747772
    ## 549 lh_cluster_2   lh      129460  -2.669875
    ## 550 lh_cluster_2   lh      129468  -2.778847
    ## 551 lh_cluster_2   lh      129981  -2.826816
    ## 552 lh_cluster_2   lh      129987  -2.769923
    ## 553 lh_cluster_2   lh      130923  -2.736296
    ## 554 lh_cluster_2   lh      130932  -2.824289
    ## 555 lh_cluster_2   lh      131061  -2.770342
    ## 556 lh_cluster_2   lh      131213  -1.856283
    ## 557 lh_cluster_2   lh      131232  -2.187672
    ## 558 lh_cluster_2   lh      131261  -1.713503
    ## 559 lh_cluster_2   lh      131862  -2.120339
    ## 560 lh_cluster_2   lh      132829  -2.178810
    ## 561 lh_cluster_2   lh      132903  -2.196242
    ## 562 lh_cluster_2   lh      132937  -2.962943
    ## 563 lh_cluster_2   lh      132947  -2.777386
    ## 564 lh_cluster_2   lh      132955  -2.883126
    ## 565 lh_cluster_2   lh      132957  -2.922586
    ## 566 lh_cluster_2   lh      132958  -3.126752
    ## 567 lh_cluster_2   lh      133048  -1.903787
    ## 568 lh_cluster_2   lh      133050  -1.896930
    ## 569 lh_cluster_2   lh      133059  -1.752161
    ## 570 lh_cluster_2   lh      133067  -1.670437
    ## 571 lh_cluster_2   lh      133266  -2.295657
    ## 572 lh_cluster_2   lh      133318  -2.077978
    ## 573 lh_cluster_2   lh      133353  -2.748770
    ## 574 lh_cluster_2   lh      133363  -2.778066
    ## 575 lh_cluster_2   lh      133388  -2.457157
    ## 576 lh_cluster_2   lh      133396  -2.388834
    ## 577 lh_cluster_2   lh      137766  -2.465033
    ## 578 lh_cluster_2   lh      137801  -1.855743
    ## 579 lh_cluster_2   lh      137853  -3.067700
    ## 580 lh_cluster_2   lh      137867  -2.974158
    ## 581 lh_cluster_2   lh      137869  -2.929076
    ## 582 lh_cluster_2   lh      137930  -1.877492
    ## 583 lh_cluster_2   lh      137942  -2.456002
    ## 584 lh_cluster_2   lh      137947  -2.299135
    ## 585 lh_cluster_2   lh      137964  -2.365802
    ## 586 lh_cluster_2   lh      138030  -2.259001
    ## 587 lh_cluster_2   lh      138050  -1.972327
    ## 588 lh_cluster_2   lh      138071  -1.909967
    ## 589 lh_cluster_2   lh      138078  -2.105213
    ## 590 lh_cluster_2   lh      138083  -1.941128
    ## 591 lh_cluster_2   lh      138090  -1.823496
    ## 592 lh_cluster_2   lh      138677  -2.668877
    ## 593 lh_cluster_2   lh      138680  -2.613598
    ## 594 lh_cluster_2   lh      138703  -2.625983
    ## 595 lh_cluster_2   lh      138712  -2.634187
    ## 596 lh_cluster_2   lh      139548  -2.139704
    ## 597 lh_cluster_2   lh      139571  -2.116695
    ## 598 lh_cluster_2   lh      139576  -2.215816
    ## 599 lh_cluster_2   lh      139680  -2.473578
    ## 600 lh_cluster_2   lh      139721  -3.084668
    ## 601 lh_cluster_2   lh      139723  -3.114415
    ## 602 lh_cluster_2   lh      139727  -3.065696
    ## 603 lh_cluster_2   lh      139738  -3.157661
    ## 604 lh_cluster_2   lh      140010  -1.873192
    ## 605 lh_cluster_2   lh      140022  -2.428306
    ## 606 lh_cluster_2   lh      140053  -2.360365
    ## 607 lh_cluster_2   lh      140081  -2.549805
    ## 608 lh_cluster_2   lh      140091  -2.649093
    ## 609 lh_cluster_2   lh      140095  -2.699740
    ## 610 lh_cluster_2   lh      140138  -2.578142
    ## 611 lh_cluster_2   lh      140152  -2.367616
    ## 612 lh_cluster_2   lh      144364  -1.877462
    ## 613 lh_cluster_2   lh      144392  -2.797173
    ## 614 lh_cluster_2   lh      144395  -2.834889
    ## 615 lh_cluster_2   lh      144409  -2.667339
    ## 616 lh_cluster_2   lh      144427  -1.869190
    ## 617 lh_cluster_2   lh      144442  -1.802146
    ## 618 lh_cluster_2   lh      144449  -2.543071
    ## 619 lh_cluster_2   lh      144451  -2.584459
    ## 620 lh_cluster_2   lh      144513  -1.768680
    ## 621 lh_cluster_2   lh      144526  -1.876261
    ## 622 lh_cluster_2   lh      144540  -1.740638
    ## 623 lh_cluster_2   lh      144651  -2.649130
    ## 624 lh_cluster_2   lh      144653  -2.673120
    ## 625 lh_cluster_2   lh      144678  -3.003641
    ## 626 lh_cluster_2   lh      144682  -3.030228
    ## 627 lh_cluster_2   lh      144684  -2.988099
    ## 628 lh_cluster_2   lh      144691  -3.092796
    ## 629 lh_cluster_2   lh      144697  -3.114132
    ## 630 lh_cluster_2   lh      144699  -3.025141
    ## 631 lh_cluster_2   lh      144705  -3.005927
    ## 632 lh_cluster_2   lh      144718  -3.012284
    ## 633 lh_cluster_2   lh      144733  -3.074111
    ## 634 lh_cluster_2   lh      144752  -3.081471
    ## 635 lh_cluster_2   lh      144771  -2.844472
    ## 636 lh_cluster_2   lh      144785  -2.831925
    ## 637 lh_cluster_2   lh      144807  -2.809669
    ## 638 lh_cluster_2   lh      144811  -2.759856
    ## 639 lh_cluster_2   lh      144917  -1.705038
    ## 640 lh_cluster_2   lh      144975  -1.935498
    ## 641 lh_cluster_2   lh      145018  -2.004066
    ## 642 lh_cluster_2   lh      145019  -1.960363
    ## 643 lh_cluster_2   lh      145074  -2.692695
    ## 644 lh_cluster_2   lh      145087  -2.799384
    ## 645 lh_cluster_2   lh      145090  -2.845287
    ## 646 lh_cluster_2   lh      145094  -2.894552
    ## 647 lh_cluster_2   lh      145188  -2.306647
    ## 648 lh_cluster_2   lh      145225  -1.855943
    ## 649 lh_cluster_2   lh      149824  -2.013308
    ## 650 lh_cluster_2   lh      150004  -1.877636
    ## 651 lh_cluster_2   lh      150009  -1.840493
    ## 652 lh_cluster_2   lh      150076  -2.418087
    ## 653 lh_cluster_2   lh      150112  -2.223941
    ## 654 lh_cluster_2   lh      150155  -3.594900
    ## 655 lh_cluster_2   lh      150161  -3.420627
    ## 656 lh_cluster_2   lh      150173  -2.852933
    ## 657 lh_cluster_2   lh      150260  -2.888609
    ## 658 lh_cluster_2   lh      150272  -2.915713
    ## 659 lh_cluster_2   lh      150275  -2.899922
    ## 660 lh_cluster_2   lh      150279  -2.885434
    ## 661 lh_cluster_2   lh      150301  -3.127359
    ## 662 lh_cluster_2   lh      150304  -3.053297
    ## 663 lh_cluster_2   lh      150310  -2.841323
    ## 664 lh_cluster_2   lh      150399  -2.940142
    ## 665 lh_cluster_2   lh      150421  -1.964484
    ## 666 lh_cluster_2   lh      150494  -3.081033
    ## 667 lh_cluster_2   lh      151055  -1.668231
    ## 668 lh_cluster_2   lh      151166  -2.402683
    ## 669 lh_cluster_2   lh      151228  -2.136434
    ## 670 lh_cluster_2   lh      151259  -2.220501
    ## 671 lh_cluster_2   lh      151311  -2.758924
    ## 672 lh_cluster_2   lh      151314  -2.776259
    ## 673 lh_cluster_2   lh      151413  -2.767905
    ## 674 lh_cluster_2   lh      155934  -2.029620
    ## 675 lh_cluster_1   lh        4821   2.963625
    ## 676 lh_cluster_1   lh        5425   1.997769
    ## 677 lh_cluster_1   lh        9117   2.490170
    ## 678 lh_cluster_1   lh        9841   2.063780
    ## 679 lh_cluster_1   lh       11366   2.144142
    ## 680 lh_cluster_1   lh       12627   2.008087
    ## 681 lh_cluster_1   lh       22057   2.075435
    ## 682 lh_cluster_1   lh       22060   1.963660
    ## 683 lh_cluster_1   lh       22065   1.781243
    ## 684 lh_cluster_1   lh       22066   1.780069
    ## 685 lh_cluster_1   lh       22067   1.807639
    ## 686 lh_cluster_1   lh       23990   1.910698
    ## 687 lh_cluster_1   lh       28193   1.899689
    ## 688 lh_cluster_1   lh       28211   2.097919
    ## 689 lh_cluster_1   lh       29200   2.083483
    ## 690 lh_cluster_1   lh       30848   2.565757
    ## 691 lh_cluster_1   lh       30851   1.985112
    ## 692 lh_cluster_1   lh       31010   1.854087
    ## 693 lh_cluster_1   lh       33265   2.018873
    ## 694 lh_cluster_1   lh       37564   1.845441
    ## 695 lh_cluster_1   lh       37582   2.767477
    ## 696 lh_cluster_1   lh       37594   2.417317
    ## 697 lh_cluster_1   lh       37601   2.274525
    ## 698 lh_cluster_1   lh       39777   2.397180
    ## 699 lh_cluster_1   lh       45741   2.006130
    ## 700 lh_cluster_1   lh       46226   2.439538
    ## 701 lh_cluster_1   lh       47540   2.152936
    ## 702 lh_cluster_1   lh       50933   1.855518
    ## 703 lh_cluster_1   lh       50938   2.877723
    ## 704 lh_cluster_1   lh       52757   2.076313
    ## 705 lh_cluster_1   lh       52781   1.998772
    ## 706 lh_cluster_1   lh       57017   1.965328
    ## 707 lh_cluster_1   lh       60568   1.815858
    ## 708 lh_cluster_1   lh       60916   2.006209
    ## 709 lh_cluster_1   lh       62357   2.056747
    ## 710 lh_cluster_1   lh       67749   1.996794
    ## 711 lh_cluster_1   lh       69551   2.376573
    ## 712 lh_cluster_1   lh       69558   2.451164
    ## 713 lh_cluster_1   lh       73522   2.716175
    ## 714 lh_cluster_1   lh       73524   2.764723
    ## 715 lh_cluster_1   lh       73528   2.547196
    ## 716 lh_cluster_1   lh       83668   1.821302
    ## 717 lh_cluster_1   lh       83740   2.489569
    ## 718 lh_cluster_1   lh       89147   2.045287
    ## 719 lh_cluster_1   lh       89203   2.460077
    ## 720 lh_cluster_1   lh       89220   1.969035
    ## 721 lh_cluster_1   lh       89223   1.948439
    ## 722 lh_cluster_1   lh       89227   1.816426
    ## 723 lh_cluster_1   lh       94116   1.932315
    ## 724 lh_cluster_1   lh       94191   2.072604
    ## 725 lh_cluster_1   lh       97316   1.963499
    ## 726 lh_cluster_1   lh       97325   1.700736
    ## 727 lh_cluster_1   lh      100763   1.929818
    ## 728 lh_cluster_1   lh      100767   2.004781
    ## 729 lh_cluster_1   lh      101638   2.424464
    ## 730 lh_cluster_1   lh      110103   1.764825
    ## 731 lh_cluster_1   lh      110142   2.678272
    ## 732 lh_cluster_1   lh      110164   2.095730
    ## 733 lh_cluster_1   lh      113392   2.018968
    ## 734 lh_cluster_1   lh      113426   2.341359
    ## 735 lh_cluster_1   lh      113437   1.979370
    ## 736 lh_cluster_1   lh      113440   1.792273
    ## 737 lh_cluster_1   lh      117414   1.957137
    ## 738 lh_cluster_1   lh      119928   2.003705
    ## 739 lh_cluster_1   lh      124674   1.854382
    ## 740 lh_cluster_1   lh      124676   1.819994
    ## 741 lh_cluster_1   lh      124693   3.034767
    ## 742 lh_cluster_1   lh      124710   2.346529
    ## 743 lh_cluster_1   lh      124713   2.297761
    ## 744 lh_cluster_1   lh      124726   2.009181
    ## 745 lh_cluster_1   lh      127229   1.925664
    ## 746 lh_cluster_1   lh      127232   1.958064
    ## 747 lh_cluster_1   lh      127234   2.055418
    ## 748 lh_cluster_1   lh      129885   1.925981
    ## 749 lh_cluster_1   lh      134634   1.817920
    ## 750 lh_cluster_1   lh      134671   1.918198
    ## 751 lh_cluster_1   lh      134674   1.970797
    ## 752 lh_cluster_1   lh      136499   1.985039
    ## 753 lh_cluster_1   lh      136504   2.051935
    ## 754 lh_cluster_1   lh      136517   2.007598
    ## 755 lh_cluster_1   lh      141347   1.929966
    ## 756 lh_cluster_1   lh      141377   2.666265
    ## 757 lh_cluster_1   lh      153697   2.753929
    ## 758 lh_cluster_1   lh      153701   2.992773
    ## 759 lh_cluster_1   lh      153710   2.540332
    ## 760 lh_cluster_1   lh      153722   2.390890
    ## 761 lh_cluster_1   lh      153727   2.651097
    ## 762 lh_cluster_1   lh      160228   2.073562
    ## 763 lh_cluster_1   lh      160232   1.979324
    ## 764 lh_cluster_1   lh      160235   2.088021
    ## 765 lh_cluster_1   lh      160302   2.003011
    ## 766 lh_cluster_1   lh      160304   2.030677
    ## 767 lh_cluster_1   lh      160322   2.036019
    ## 768 lh_cluster_1   lh      160334   1.702631
    ## 769 lh_cluster_1   lh      160338   1.757455

This concludes the package vignette, thanks for reading.
