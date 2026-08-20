# Transform MNI305 coords (FreeSurfer fsaverage surface) to MNI152 coordinates.

Transform MNI305 coords (FreeSurfer fsaverage surface) to MNI152
coordinates.

## Usage

``` r
coord_MNI305_to_MNI152(
  vertex_coords,
  method = getOption("brainloc.method_MNI305_to_from_MNI152", default = "linear"),
  surface = "orig",
  fs_home = getOption("brainloc.fs_home", default = Sys.getenv("FREESURFER_HOME"))
)
```

## Arguments

- vertex_coords:

  `nx3` numerical matrix of MNI305 surface RAS coordinates, typically
  from fsaverage surface vertices.

- method:

  character string, the method to use to map from MNI305 to MNI152 along
  the way. One of "best_available", "regfusionr", and "linear".

- surface:

  optional character string or
  [`hemilist`](https://dfsp-spirit.github.io/brainloc/reference/hemilist.md)
  of surfaces, the surface to use to find a surface vertex close to the
  given query coordinates. Only used if 'method' results in `regfusionr`
  being used. Passed on to
  [`regfusionr::mni305_coords_to_mni152_coords`](https://rdrr.io/pkg/regfusionr/man/mni305_coords_to_mni152_coords.html).

- fs_home:

  optional character string, the path to the FREESURFER_HOME directory
  from which to load the surfaces from the 'surface' parameter. Only
  used if 'method' results in `regfusionr` being used. Passed on to
  [`regfusionr::mni305_coords_to_mni152_coords`](https://rdrr.io/pkg/regfusionr/man/mni305_coords_to_mni152_coords.html).

## Value

`nx3` numerical matrix of MNI152 coords.

## Note

To verify MNI305→MNI152 results, you can click a vertex in FreeView to
get its MNI305 coordinates, run them through this function, and compare
the MNI152 output with the values shown by the BioImageSuite
MNI→Talairach tool at
`https://bioimagesuiteweb.github.io/bisweb-manual/tools/mni2tal.html`.
Note that this external tool uses a different Talairach transform
(Lacadie et al.) than the one used by
[`coord_MNI152_to_talairach`](https://dfsp-spirit.github.io/brainloc/reference/coord_MNI152_to_talairach.md)
(Brett), so MNI152 values should match but Talairach values may differ
slightly.

## Examples

``` r
if (FALSE) { # \dontrun{
# Get MNI152 coordinates for first 3 fsaverage lh vertices:
surf = freesurferformats::read.fs.surface("/opt/freesurfer/subjects/fsaverage/surf/lh.white");
coord_MNI305_to_MNI152(surf$vertices[1:3, ]);
} # }
```
