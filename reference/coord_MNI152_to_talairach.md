# Transform MNI152 coords to approximate Talairach space using Matthew Brett's approach.

Transform MNI152 coordinates to Talairach space using the piecewise
linear transform described by Matthew Brett at
`http://brainmap.org/training/BrettTransform.html`. This is an
\*\*approximation\*\* of Talairach space, not the original Talairach
atlas space (which is defined by the AC-PC line and a bounding box from
a single post-mortem brain). The transform applies different scaling
factors above and below the AC line to account for differences between
the MNI152 template brain and the Talairach atlas brain.

## Usage

``` r
coord_MNI152_to_talairach(mni152_coords)
```

## Arguments

- mni152_coords:

  nx3 numerical matrix of RAS coordinates in MNI152 space.

## Value

nx3 numerical matrix of approximate Talairach coordinates.

## Note

When using the output of this function with
[`get_talairach_label`](https://dfsp-spirit.github.io/brainloc/reference/get_talairach_label.md)
to query the Talairach Daemon atlas (`talairach.org`), be aware that the
Daemon is defined in true Talairach space while this transform produces
approximate Talairach coordinates. In practice, the Talairach Daemon has
fairly coarse labels (lobe-level), so small coordinate discrepancies
usually do not affect the result. This approach is widely used in the
field (SPM, FSL, etc.).

This implementation is published under the GPL license. See
`https://github.com/sccn/dipfit/blob/master/mni2tal_matrix.m` and
`https://github.com/sccn/dipfit/blob/master/mni2tal.m` for a Matlab
implementation of the method. All credits go to Matthew Brett.

## Examples

``` r
    mni_coords = matrix(c(10, 12, 14), nrow = 1, ncol = 3, byrow = TRUE);
    coord_MNI152_to_talairach(mni_coords);
#>      [,1]     [,2]     [,3]
#> [1,]  9.9 12.27003 12.28254
```
