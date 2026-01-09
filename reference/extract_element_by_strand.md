# Extract Genomic Elements by Strand

Filters a data frame or GTF/GFF3 file to extract specific genomic
elements (genes, transcripts, exons, or introns) located on a specified
DNA strand (+ or -).

## Usage

``` r
extract_element_by_strand(input, type, strand, verbose = TRUE)
```

## Arguments

- input:

  A data frame derived from a GTF/GFF3 file or containing genomic
  annotations.

- type:

  A character string specifying the type of element to extract. Valid
  options: `"gene"`, `"transcript"`, `"exon"`, `"intron"`.

- strand:

  A character string indicating the DNA strand to filter. Valid options:
  `"+"` (forward) or `"-"` (reverse).

- verbose:

  A logical indicating whether to print the count of extracted elements.
  Defaults to `TRUE`.

## Value

A data frame containing only the elements of the specified type located
on the chosen strand.

## Details

This function filters the input data based on the `type` and `strand`
columns. It is useful for strand-specific analyses, such as studying
antisense transcripts or strand-biased genomic features.

## Examples

``` r
# Example 1: Extract genes on the forward strand using gencode.v1
file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
gtf_v1 <- load_file(file_v1)
forward_genes_v1 <- extract_element_by_strand(gtf_v1, type = "gene", strand = "+")
#> Genes on strand '+': 3

# Example 2: Extract exons on the reverse strand using gencode.v1
reverse_exons_v1 <- extract_element_by_strand(gtf_v1, type = "exon", strand = "-")
#> Exons on strand '-': 8
```
