# Create a brain parcellation from a subject in a FreeSurfer recon-all output directory

Create a brain parcellation from a subject in a FreeSurfer recon-all
output directory

## Usage

``` r
brainparc_fs(subjects_dir, subject_id, surface = "white", atlas = c("aparc"))
```

## Arguments

- subjects_dir:

  string. The FreeSurfer SUBJECTS_DIR, i.e., a directory containing the
  data for all your subjects, each in a subdir named after the subject
  identifier.

- subject_id:

  string. The subject identifier

- surface:

  string. The surface name. E.g., "white", or "pial". Used to construct
  the name of the surface file to be loaded.

- atlas:

  vector of character strings, the brain atlases to use. For this
  FreeSurfer version, something like "aparc" for the Desikan-Killiani
  atlas, "aparc.a2009s" for Destrieux, etc.

## Value

a `brainparc` instance
