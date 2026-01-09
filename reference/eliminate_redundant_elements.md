# Eliminate Redundant Genomic Elements

Removes redundant genomic elements (exons or introns) from a data frame,
ensuring each element is uniquely represented. Redundancy is determined
by genomic coordinates and gene ID.

## Usage

``` r
eliminate_redundant_elements(input, element_type = "exon")
```

## Arguments

- input:

  A data frame containing genomic element coordinates (exons or introns)
  with columns `seqnames`, `start`, `end`, and `gene_id`.

- element_type:

  The type of genomic element to process. Valid options are `"exon"`
  (default) or `"intron"`.

## Value

A data frame with redundant elements removed, retaining only unique
entries based on coordinates and gene ID.

## Details

This function uses genomic coordinates (chromosome, start, end) and gene
ID to identify and remove duplicate entries. For exons, coordinates are
directly compared. For introns, coordinates are derived from
`intron_start` and `intron_end` columns (check `extract_introns`
function for more details)

## Examples

``` r
file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
gtf_v1 <- load_file(file_v1)
# Eliminate redundant exons
nonredundant_exons <- eliminate_redundant_elements(gtf_v1, element_type = "exon")
```
