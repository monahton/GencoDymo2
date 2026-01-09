# Calculate GC Content of Genomic Features

Computes the GC content percentage for exons or introns using genome
sequence data from a BSgenome object.

## Usage

``` r
calculate_gc_content(input, genome, verbose = TRUE)
```

## Arguments

- input:

  A data frame containing genomic features (exons/introns) with
  `seqnames`, `start` or `intron_start`, `end` or `intron_end`, and
  `strand` columns.

- genome:

  A BSgenome object representing the reference genome (e.g.,
  `BSgenome.Hsapiens.UCSC.hg38`).

- verbose:

  A logical indicating whether to print progress messages. Defaults to
  `TRUE`.

## Value

The input data frame with an additional `gc_content` column containing
GC percentages for each feature.

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
  gtf_v1 <- load_file(file_v1)
  gtf_with_gc <- calculate_gc_content(gtf_v1, genome = genome, verbose = FALSE)
}
} # }
```
