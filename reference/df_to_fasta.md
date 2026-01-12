# Convert Data Frame to FASTA File

Converts a data frame containing sequence IDs and sequences into a
FASTA-formatted file, optionally compressed as gzip.

## Usage

``` r
df_to_fasta(df, id_col, seq_col, output_file = NULL, gzip = TRUE, verbose = TRUE)
```

## Arguments

- df:

  A data frame with at least two columns: one for sequence IDs and one
  for sequences.

- id_col:

  A character string specifying the column name containing sequence IDs.

- seq_col:

  A character string specifying the column name containing sequence
  data.

- output_file:

  A character string specifying the output file path. If NULL, the
  function will stop with an informative message.

- gzip:

  A logical indicating whether to compress the output as a gzip file.
  Defaults to `TRUE`.

- verbose:

  A logical indicating whether to print progress messages. Defaults to
  `TRUE`.

## Value

No return value. Writes a FASTA file to the specified path.

## Details

This function efficiently writes large sequence datasets to FASTA
format, handling compression and progress reporting. It validates input
columns and manages memory by processing data in chunks.

## Examples

``` r
temp_dir <- tempdir()
temp_output <- file.path(temp_dir, "output.fa.gz")
seq_data <- data.frame(
  transcript_id = c("ENST0001", "ENST0002"),
  sequence = c("ATGCTAGCTAG", "GCTAGCTAGCT")
)
df_to_fasta(seq_data, "transcript_id", "sequence", temp_output)
#> Writing FASTA file...
#> FASTA file successfully saved to /tmp/RtmpYUnJ5S/output.fa.gz 
```
