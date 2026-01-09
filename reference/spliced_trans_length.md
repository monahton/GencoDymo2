# Calculate Spliced Transcript Lengths

Computes the spliced length of transcripts by summing the lengths of
their constituent exons. This represents the mature RNA length after
intron removal.

## Usage

``` r
spliced_trans_length(input)
```

## Arguments

- input:

  A data frame containing GENCODE annotations in GTF format, including
  `type`, `transcript_id`, and `width` columns.

## Value

A data frame with two columns: `transcript_id` and `transcript_length`,
where the latter is the sum of exon widths for each transcript.

## Details

This function processes the input data to:

1.  Filter entries to include only exons.

2.  Group exons by transcript ID.

3.  Sum the widths of all exons per transcript.

The result provides the total length of the mature spliced RNA for each
transcript.

## Examples

``` r
file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
gtf_v1 <- load_file(file_v1)
spliced_lengths <- spliced_trans_length(gtf_v1)
```
