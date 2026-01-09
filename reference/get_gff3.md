# Download GFF3 File from the GENCODE Database

Downloads a GFF3 file for a specified species, release version, and
annotation type from the GENCODE database. The file is saved to a
user-specified directory or the current working directory by default.

## Usage

``` r
get_gff3(species, release_version, annotation_type, dest_folder)
```

## Arguments

- species:

  A character string indicating the species. Supported values are:

  - "human"

  - "mouse"

- release_version:

  A character string specifying the release version. Options include:

  - "latest_release": Fetches the latest release for the species.

  - "release_X": Specific release version for human (e.g.,
    "release_42").

  - "release_MX": Specific release version for mouse (e.g.,
    "release_M36").

- annotation_type:

  A character string specifying the annotation type. Supported values
  include:

  - "annotation.gff3.gz"

  - "basic.annotation.gff3.gz"

  - "chr_patch_hapl_scaff.annotation.gff3.gz"

  - "chr_patch_hapl_scaff.basic.annotation.gff3.gz"

  - "long_noncoding_RNAs.gff3.gz"

  - "primary_assembly.annotation.gff3.gz"

  - "primary_assembly.basic.annotation.gff3.gz"

  - "tRNAs.gff3.gz"

  - "polyAs.gff3.gz"

- dest_folder:

  A character string specifying the destination folder. Defaults to the
  current working directory.

## Value

A character string specifying the full path of the downloaded GFF3 file.

## Details

The function dynamically determines the correct file URL based on the
provided parameters and downloads the GFF3 file to the desired location.
If "latest_release" is specified for `release_version`, the function
will first determine the latest available release using
[`get_latest_release()`](https://monahton.github.io/GencoDymo2/reference/get_latest_release.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Download the latest human GTF file with primary assembly annotations into a temp directory
temp_dir <- tempdir()
gff3_file <- get_gff3(
  species = "human",
  release_version = "latest_release",
  annotation_type = "primary_assembly.basic.annotation.gff3.gz",
  dest_folder = temp_dir
)
print(gff3_file)

# Download a specific mouse release with long noncoding RNA annotations into a temp directory
temp_dir <- tempdir()
gff3_file <- get_gff3(
  species = "mouse",
  release_version = "release_M36",
  annotation_type = "long_noncoding_RNAs.gff3.gz",
  dest_folder = temp_dir
)
print(gff3_file)
} # }
```
