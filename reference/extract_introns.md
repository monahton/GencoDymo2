# Extract Intron Coordinates from GENCODE Annotations

Processes a GTF file or data frame to extract intron coordinates,
including their genomic positions, transcript associations, and
metadata. The function handles both positive and negative strands,
ensuring correct orientation of intron boundaries. It is designed to
work with GENCODE-formatted annotations.

## Usage

``` r
extract_introns(input, verbose = TRUE)
```

## Arguments

- input:

  A character string specifying the file path to a GTF/GFF3 file, or a
  data frame containing GTF data. The input must include columns:
  `seqnames`, `start`, `end`, `strand`, `type`, and `transcript_id`.

- verbose:

  Logical. If `TRUE` (default), progress messages are printed during
  execution.

## Value

A data frame with the following columns:

- `seqnames`: Chromosome or scaffold name.

- `intron_start`, `intron_end`: Genomic start/end positions of the
  intron.

- `width`: Length of the intron (intron_end - intron_start + 1).

- `gene_id`, `transcript_id`, `intron_number`: Gene/transcript
  identifiers and intron position within the transcript.

- Additional metadata columns from the original GTF (e.g., `gene_name`).

## Details

1.  Input Handling: Loads the GTF file if a path is provided; uses the
    data frame directly if supplied.

2.  Exon Filtering: Extracts exon records and classifies them (e.g.,
    first/last exon) if needed.

3.  Multi-Exon Transcripts: Removes single-exon transcripts to focus on
    spliced transcripts.

4.  Intron Calculation: Determines intron coordinates by identifying
    gaps between consecutive exons, adjusting for strand orientation.

5.  Output: Returns a data frame with intron coordinates and metadata,
    sorted by gene and position.

## See also

[`load_file`](https://monahton.github.io/GencoDymo2/reference/load_file.md),
[`classify_exons`](https://monahton.github.io/GencoDymo2/reference/classify_exons.md)

## Examples

``` r
file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")

# From a pre-loaded data frame
gtf_v1 <- load_file(file_v1)
introns <- extract_introns(gtf_v1)
#> Preparing input...
#> Removing single-exon transcripts...
#> Single-exon transcripts: 1
#> Extracting intron coordinates...
#> Collecting intron data...
#> Total introns: 3
```
