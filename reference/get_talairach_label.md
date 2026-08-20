# Retrieve Talairach labels for Talairach coordinates from `talairach.org` files.

Retrieve Talairach labels for Talairach coordinates from `talairach.org`
files.

## Usage

``` r
get_talairach_label(
  tal_coords,
  talairach_vol_file = NULL,
  lookup_table_file = NULL,
  check_oob = TRUE
)
```

## Arguments

- tal_coords:

  `nx3` numeric matrix, the `n` query Talairach coordinates for which
  you want to retrieve the labels. (Three columns, each row represents
  the x,y,z coordinates of a point in Talairach space.)

- talairach_vol_file:

  optional character string, the path to the NIFTI Talairach label
  volume file from `talairach.org`. It is named `talairach.nii` on the
  website. You can also gzip it and use a much smaller `.nii.gz` version
  of the file. Leave at `NULL` to auto-download from `talairach.org` if
  needed.

- lookup_table_file:

  optional character string, the path to the Talairach label index text
  file from `talairach.org`. It is named `labels.txt` on the website.
  Leave at `NULL` to auto-download from `talairach.org` if needed.

- check_oob:

  logical, whether to check for out-of-bounds voxels, which do not fall
  into the NIFTI volume during the computation (because the talairach
  coordinate is a bit off). Checking is slower, but if such coordinates
  exist in your query and you do not check, the function will stop with
  an error. Leave this alone if in doubt.

## Value

data.frame describing labels for the coordinates. The following columns
are included: `cx,cy,cz`: the query coordinates, as given in parameter
'tal_coords'. `v1,v2,v3`: The voxel indices in the talairach.nii file
that map to the query coordinates. `label_lvl1,...,label_lvl5`: the
label strings for the location, in a hierarchy with 5 levels (parsed
from `label_full` for you as a convenience). `label_full`: The raw, full
label string for the location. A star `*` means that no label name is
available for the given coordinate at this level.

## Note

If you leave the file path arguments `talairach_vol_file` and
`lookup_table_file` at NULL, the files will be downloaded from
`talairach.org` automatically if that is possible (e.g., you have a
working internet connection and the server is up). If it fails, the
function will stop in the next step.

When using this function or the Talairach.org data in your research,
please cite the following two publications:
`Lancaster JL, Woldorff MG, Parsons LM, Liotti M, Freitas CS, Rainey L, Kochunov PV, Nickerson D, Mikiten SA, Fox PT, "Automated Talairach Atlas labels for functional brain mapping". Human Brain Mapping 10:120-131, 2000.`
and
`Lancaster JL, Rainey LH, Summerlin JL, Freitas CS, Fox PT, Evans AC, Toga AW, Mazziotta JC. Automated labeling of the human brain: A preliminary report on the development and evaluation of a forward-transform method. Hum Brain Mapp 5, 238-242, 1997.`

Currently it is not possible to search 'in the vicinity' to compensate
for a slightly inaccurate query coordinate that may be outside of the
brain by just a tiny bit. Please open an issue if you feel this is
required.

## Examples

``` r
if (FALSE) { # \dontrun{
query_tal_coords = matrix(seq.int(15), nrow=5, ncol=3);
get_talairach_label(query_tal_coords);
} # }
```
