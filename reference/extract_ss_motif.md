# Extract Splice Site Motifs for MaxEntScan Analysis (5' or 3')

Extracts donor (5') or acceptor (3') splice site motifs from a genomic
dataset using BSgenome sequences.

## Usage

``` r
extract_ss_motif(input, genome, type, verbose, save_fasta, output_file)
```

## Arguments

- input:

  A data frame with columns: `seqnames`, `strand`, `intron_start`,
  `intron_end`, `transcript_id`, `intron_number`.

- genome:

  A BSgenome object (e.g., from BSgenome.Hsapiens.UCSC.hg38).

- type:

  One of `"5ss"` (donor) or `"3ss"` (acceptor).

- verbose:

  Logical; print progress messages.

- save_fasta:

  Logical; write a FASTA file of extracted motifs.

- output_file:

  Name/path of the output FASTA file.

## Value

A data frame including extracted splice site motifs.

## See also

[`assign_splice_sites`](https://monahton.github.io/GencoDymo2/reference/assign_splice_sites.md),
[`df_to_fasta`](https://monahton.github.io/GencoDymo2/reference/df_to_fasta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
  gtf_v1 <- load_file(file_v1)
  introns <- extract_introns(gtf_v1)
  motifs_df <- extract_ss_motif(introns, genome, "5ss", verbose = FALSE)
}
} # }
```
