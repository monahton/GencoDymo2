# Classify Exons by Their Relative Transcript Position

Categorizes each exon in a transcript as single, first, inner, or last
based on its position. This classification helps in analyzing transcript
structure and splicing patterns.

## Usage

``` r
classify_exons(input, verbose = TRUE)
```

## Arguments

- input:

  A character string specifying the path to a GTF/GFF3 file or a data
  frame containing GTF-formatted data. The file should follow the
  GENCODE format conventions.

- verbose:

  A logical indicating whether to display progress messages and exon
  counts. Defaults to `TRUE`.

## Value

A data frame identical to the input but with an additional column
`EXON_CLASSIFICATION`. This column contains one of the following values
for each exon:

- `"single_exon"`: The transcript contains only one exon.

- `"first_exon"`: The first exon in a multi-exon transcript.

- `"last_exon"`: The last exon in a multi-exon transcript.

- `"inner_exon"`: Exons between the first and last in multi-exon
  transcripts.

## Details

The function processes the input GTF data to:

1.  Load the GTF file (if a path is provided) or use the provided data
    frame.

2.  Filter entries to include only exons.

3.  For each transcript, count the total exons and determine each exon's
    position.

4.  Classify exons based on their position and the total count.

If `verbose = TRUE`, the function prints counts of each exon type.

## Examples

``` r
# Example 1: Using the provided example GTF files
file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
gtf_v1 <- classify_exons(file_v1)
#> Loading GTF file: /home/runner/work/_temp/Library/GencoDymo2/extdata/gencode.v1.example.gtf.gz
#> Classifying exons...
#> Single exons: 1
#> First exons: 2
#> Last exons: 2
#> Inner exons: 1

# Example 2: Using a pre-loaded data frame
file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
gtf_v1 <- load_file(file_v1)
classified_data_v1 <- classify_exons(gtf_v1)
#> Using input dataframe...
#> Classifying exons...
#> Single exons: 1
#> First exons: 2
#> Last exons: 2
#> Inner exons: 1
```
