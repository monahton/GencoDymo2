# Generate Summary Statistics for Genomic Elements

Calculates descriptive summary statistics (mean, median, standard
deviation, etc.) for the lengths of exons or introns, grouped by
classification (e.g., first exon, inner intron).

## Usage

``` r
stat_summary(input, type, verbose = TRUE)
```

## Arguments

- input:

  A data frame containing classified exons or introns. For exons, must
  include `EXON_CLASSIFICATION`. For introns, requires classification
  columns.

- type:

  A character string specifying the element type. Valid options:
  `"exon"` or `"intron"`.

- verbose:

  A logical indicating whether to print progress messages. Defaults to
  `TRUE`.

## Value

A data frame with summary statistics for each element group, including
mean, median, standard deviation, standard error, quartiles, and sample
size.

## Details

For exons, statistics are grouped by `EXON_CLASSIFICATION`. For introns,
groups include `first_intron`, `inner_intron`, and splice site types
(e.g., `gc_intron`).

## Examples

``` r
file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
gtf_v1 <- load_file(file_v1)
# Exon statistics
exon_stats <- stat_summary(gtf_v1, type = "exon")
#> Classifying exons...
#> Exon classification completed.
#> Calculating summary statistics...
#> Summary statistics calculation completed.

```
