# Identify Single-Exon Genes/Transcripts in GENCODE Data

This function identifies single-exon genes or transcripts from a genomic
annotation dataset, typically obtained from the GENCODE database. It
processes input data, either from a file or an already loaded data
frame, and returns either single-exon genes or single-exon transcripts
based on the user's selection. The function works seamlessly with
GTF/GFF files or data frames that include the standard GENCODE
annotation fields. This is especially useful for identifying genes or
transcripts that do not undergo splicing.

## Usage

``` r
extract_single_exon(input, level = "gene", output_file = NULL)
```

## Arguments

- input:

  Either a file path to a GTF/GFF3 file or a data frame representing
  genomic annotations. The input data should contain at least the
  columns `type`, `gene_id`, `transcript_id`, and `exon_number`. If a
  file path is provided, it is assumed to be a GTF/GFF file.

- level:

  A character string specifying the level of analysis. It should be
  either `"gene"` or `"transcript"`. The default is "gene". The
  selection allows users to specify whether they want to identify
  single-exon genes or single-exon transcripts in the annotation.

- output_file:

  Optional. A file path to save the results as a tab-separated file. If
  not provided, the results will not be saved. The file will include the
  IDs of the single-exon genes or transcripts along with their
  respective exon IDs. This feature is useful for exporting results for
  further analysis or reporting purposes.

## Value

A data frame containing the IDs of single-exon genes or transcripts
based on the selected level. The returned data frame contains the
following columns:

- `gene_id` or `transcript_id`: The IDs of the single-exon genes or
  transcripts.

- `exon_count`: The count of exons for each entity (always 1 for
  single-exon entities).

- `exon_id`: The ID of the exon associated with the single-exon gene or
  transcript.

The data frame can then be used for downstream analysis, such as
identifying unique single-exon genes or transcripts, or for exporting to
a file.

## Details

1.  Input Validation: Checks for required columns and valid `level`
    specification.

2.  Exon Counting: Groups exons by gene/transcript and counts
    occurrences.

3.  Single-Exon Filtering: Retains only entities with exactly one exon.

4.  Output: Returns filtered results; optionally writes to file.

## See also

[`extract_introns`](https://monahton.github.io/GencoDymo2/reference/extract_introns.md),
[`load_file`](https://monahton.github.io/GencoDymo2/reference/load_file.md)

## Examples

``` r
# Example input data frame
input <- data.frame(
  type = c("exon", "exon", "exon", "exon", "exon"),
  gene_id = c("gene1", "gene1", "gene2", "gene3", "gene4"),
  transcript_id = c("tx1", "tx1", "tx2", "tx3", "tx4"),
  exon_number = c(1, 2, 1, 1, 1),
  exon_id = c("e1", "e2", "e1", "e1", "e1")
)

# Identify single-exon genes
single_exon_genes <- extract_single_exon(input, level = "gene")
print(single_exon_genes)
#>   gene_id exon_count exon_id
#> 1   gene2          1      e1
#> 2   gene3          1      e1
#> 3   gene4          1      e1

# Identify single-exon transcripts
single_exon_transcripts <- extract_single_exon(input, level = "transcript")
print(single_exon_transcripts)
#>   transcript_id exon_count exon_id
#> 1           tx2          1      e1
#> 2           tx3          1      e1
#> 3           tx4          1      e1
```
