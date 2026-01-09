# Assign intron donor and acceptor splice sites consensus

This function takes a data frame of intron coordinates and a genome
sequence (ideally human or mouse) and returns a data frame with two
additional columns for the donor and acceptor splice site consensus
sequences. It prepares the donor and acceptor sequences based on the
provided intron coordinates and the specified genome (e.g., human hg38),
making it useful for downstream analysis of splicing events.

## Usage

``` r
assign_splice_sites(input, genome, verbose = TRUE)
```

## Arguments

- input:

  A data frame containing intron coordinates.

- genome:

  A `BSgenome` object like `BSgenome.Hsapiens.UCSC.hg38`. Must be
  explicitly passed.

- verbose:

  Logical. If TRUE, the function prints progress messages while
  preparing the splice site data. Default is TRUE.

## Value

A data frame containing the original intron data, with two additional
columns:

- `donor_ss`: The donor splice site consensus sequence for each intron.

- `acceptor_ss`: The acceptor splice site consensus sequence for each
  intron.

## See also

[`extract_introns`](https://monahton.github.io/GencoDymo2/reference/extract_introns.md),
[`find_cryptic_splice_sites`](https://monahton.github.io/GencoDymo2/reference/find_cryptic_splice_sites.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
    gtf_v1 <- load_file(file_v1)
    introns_df <- extract_introns(gtf_v1)
    result <- assign_splice_sites(introns_df, genome)
  }
} # }
```
