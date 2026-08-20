# Create a brain parcellation from custom data.

Create a brain parcellation from custom data.

## Usage

``` r
brainparc(surfaces, annots)
```

## Arguments

- surfaces:

  a named list with a single entry, which must contain a hemilist of
  fs.surface instances. The name can be anything describing the surface.
  Typically it is something like 'white', 'pial', or 'inflated' for
  FreeSurfer surfaces.

- annots:

  a named list with at least one entry, all of which must contain a
  hemilist of minannot (or `fs.annot`) instances. The name can be
  anything describing the atlas from which the parcellation originates.
  Typically it is something like 'aparc', 'aparc.a2009s', or 'Desikan'
  for FreeSurfer parcellations.

## Value

a `brainparc` instance

## Note

A `minannot` for a surface is just a vector with length equal to the
number of surface vertices, that holds one area label per vertex. An
area label is an arbitrary character string, but typically a brain
region name.
