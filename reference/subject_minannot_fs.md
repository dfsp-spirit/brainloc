# Return simplified annotation: the region names only.

Return simplified annotation: the region names only.

## Usage

``` r
subject_minannot_fs(subjects_dir, subject_id, atlas)
```

## Arguments

- subjects_dir:

  string. The FreeSurfer SUBJECTS_DIR, i.e., a directory containing the
  data for all your subjects, each in a subdir named after the subject
  identifier.

- subject_id:

  string. The subject identifier

- atlas:

  string. The atlas name. E.g., "aparc", "aparc.2009s", or
  "aparc.DKTatlas". Used to construct the name of the annotation file to
  be loaded.
