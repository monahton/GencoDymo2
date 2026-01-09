# Compare Annotation Counts Between Two GENCODE Releases

Compares the number of annotated genomic elements (genes, transcripts,
exons, introns) between two specified GENCODE releases. The function
accepts input data as either file paths (in GTF/GFF format) or
pre-loaded data frames. It computes the absolute difference (delta), the
percentage change relative to a chosen baseline, and determines the
direction of change (increase or decrease).

## Usage

``` r
compare_release(input1, input2, type, gene_biotype = NULL, baseline = "count2")
```

## Arguments

- input1:

  A character string specifying the file path to a GTF/GFF file from the
  first GENCODE release, or a data frame containing the annotation data.

- input2:

  A character string specifying the file path to a GTF/GFF file from the
  second GENCODE release, or a data frame containing the annotation
  data.

- type:

  A character string indicating the type of genomic element to compare.
  Valid options are `"gene"`, `"transcript"`, `"exon"`, or `"intron"`.

- gene_biotype:

  An optional character string specifying a particular gene biotype to
  filter comparisons (e.g., `"protein_coding"`, `"lncRNA"`). If `NULL`
  (default), all gene types are included.

- baseline:

  A character string defining the baseline for calculating percentage
  change. Options include:

  - `"count1"`: Uses the count from the first input (release) as the
    baseline.

  - `"count2"`: Uses the count from the second input (release) as the
    baseline (default).

  - `"average"`: Uses the average of the counts from both inputs as the
    baseline.

## Value

A list with the following elements:

- `delta`: The absolute difference in the number of annotations.

- `percentage`: The percentage change relative to the selected baseline.

- `direction`: A string indicating the direction of the change
  ("increase", "decrease", or "no change").

## Details

This function processes two GENCODE releases to compare annotation
counts for a specified genomic element type. Key steps include:

1.  Input Handling: If inputs are file paths, they are loaded into data
    frames using the `load_file` function. Data frames are used
    directly.

2.  Element Filtering: If `gene_biotype` is specified, annotations are
    filtered to include only that gene biotype.

3.  Count Calculation: The number of elements (genes, transcripts, etc.)
    of the specified type is counted in each release.

4.  Delta and Percentage: The absolute difference (delta) and percentage
    change are calculated based on the chosen baseline.

5.  Direction Determination: The direction of change is determined by
    comparing counts between the two releases.

The function provides both numerical results and a formatted console
output highlighting key metrics.

## See also

[`load_file`](https://monahton.github.io/GencoDymo2/reference/load_file.md),
[`get_gtf`](https://monahton.github.io/GencoDymo2/reference/get_gtf.md),
[`get_gff3`](https://monahton.github.io/GencoDymo2/reference/get_gff3.md).

## Examples

``` r
file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
file_v2 <- system.file("extdata", "gencode.v2.example.gtf.gz", package = "GencoDymo2")
# Example 1: Using data frames with the provided example GTF files
gtf_v1 <- load_file(file_v1)
gtf_v2 <- load_file(file_v2)
comparison <- compare_release(gtf_v1, gtf_v2, type = "gene")
#> Delta: 1
#> Percentage: 33.333%
#> Change Direction: increase

```
