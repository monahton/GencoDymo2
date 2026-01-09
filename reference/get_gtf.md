# Download GTF File from the GENCODE Database

Downloads a GTF file for a specified species, release version, and
annotation type from the GENCODE database. The file is saved to a
user-specified directory or the current working directory by default.

## Usage

``` r
get_gtf(species, release_version, annotation_type, dest_folder)
```

## Arguments

- species:

  A character string indicating the species. Supported values are:

  - "human"

  - "mouse"

- release_version:

  A character string specifying the release version. Options include:

  - "latest_release": Automatically fetches the latest release for the
    specified species.

  - "release_X": Specific human release (e.g., "release_47").

  - "release_MX": Specific mouse release (e.g., "release_M36").

- annotation_type:

  A character string specifying the annotation type. Valid options are:

  - "annotation.gtf.gz"

  - "basic.annotation.gtf.gz"

  - "chr_patch_hapl_scaff.annotation.gtf.gz"

  - "chr_patch_hapl_scaff.basic.annotation.gtf.gz"

  - "long_noncoding_RNAs.gtf.gz"

  - "primary_assembly.annotation.gtf.gz"

  - "primary_assembly.basic.annotation.gtf.gz"

  - "tRNAs.gtf.gz"

  - "polyAs.gtf.gz"

- dest_folder:

  A character string specifying the destination folder where the file
  will be downloaded. Defaults to the current working directory.

## Value

A character string specifying the path to the downloaded GTF file.

## Details

The function dynamically determines the correct file URL based on the
provided parameters and downloads the GTF file to the desired location.
If "latest_release" is specified for `release_version`, the function
will first determine the latest available release using
[`get_latest_release()`](https://monahton.github.io/GencoDymo2/reference/get_latest_release.md).

## Examples

``` r
if (FALSE) { # \dontrun{
# Download the latest human GTF file with primary assembly annotations into a temp directory
temp_dir <- tempdir()
gtf_file <- get_gtf(
  species = "human",
  release_version = "latest_release",
  annotation_type = "primary_assembly.basic.annotation.gtf.gz",
  dest_folder = temp_dir
)
print(gtf_file)

# Download a specific mouse release with long noncoding RNA annotations into a temp directory
temp_dir <- tempdir()
gtf_file <- get_gtf(
  species = "mouse",
  release_version = "release_M36",
  annotation_type = "long_noncoding_RNAs.gtf.gz",
  dest_folder = temp_dir
)
print(gtf_file)
} # }
```
