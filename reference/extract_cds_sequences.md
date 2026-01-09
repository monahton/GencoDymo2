# Extract Coding Sequences (CDS) from GTF Annotations

Extracts CDS regions from a GTF annotation file or data frame using
genomic coordinates and retrieves corresponding DNA sequences from a
BSgenome reference.

## Usage

``` r
extract_cds_sequences(input, genome, save_fasta, output_file, verbose)
```

## Arguments

- input:

  A character string (GTF file path) or GRanges object containing CDS
  annotations.

- genome:

  A BSgenome object for the relevant genome (e.g.,
  BSgenome.Hsapiens.UCSC.hg38).

- save_fasta:

  Logical indicating whether to save sequences to a FASTA file.

- output_file:

  Character string specifying the FASTA output path.

- verbose:

  Logical indicating whether to print progress messages.

## Value

A data frame containing CDS annotations with corresponding sequences. If
`save_fasta = TRUE`, also writes a FASTA file.

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  suppressPackageStartupMessages(library(GenomicRanges))
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
  gtf_v1 <- load_file(file_v1) # Should return GRanges
  cds_seqs <- extract_cds_sequences(gtf_v1, genome, save_fasta = FALSE)
}
} # }
```
