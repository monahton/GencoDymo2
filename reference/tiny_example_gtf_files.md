# Tiny example GTF files

These are minimal GTF files used for examples and testing within the
package.

- **gencode.v1.example.gtf.gz** contains two genes:

  - GeneA: a single-exon, unspliced gene.

  - GeneB: a spliced gene with two transcripts and multiple exons.

- **gencode.v2.example.gtf.gz** contains the same two genes as in
  `gencode.v1.example.gtf.gz` plus:

  - GeneC: a new spliced gene with multiple transcripts and many exons.

These files are stored in the `inst/extdata/` directory and can be
accessed using
`system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")`
or
`system.file("extdata", "gencode.v2.example.gtf.gz", package = "GencoDymo2")`.

## Format

Two external GTF files.

## Source

Generated manually for testing purposes.

## See also

[`load_file`](https://monahton.github.io/GencoDymo2/reference/load_file.md)

## Examples

``` r
tiny_v1_path <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
tiny_v2_path <- system.file("extdata", "gencode.v2.example.gtf.gz", package = "GencoDymo2")

gtf1 <- load_file(tiny_v1_path)
head(gtf1)
#>   seqnames start   end width strand source       type score phase gene_id
#> 1     chr1  1000  5000  4001      + HAVANA       gene    NA    NA   GeneA
#> 2     chr1  1000  5000  4001      + HAVANA transcript    NA    NA   GeneA
#> 3     chr1  1000  5000  4001      + HAVANA       exon    NA    NA   GeneA
#> 4     chr1  6000 12000  6001      - HAVANA       gene    NA    NA   GeneB
#> 5     chr1  6000 12000  6001      - HAVANA transcript    NA    NA   GeneB
#> 6     chr1  6000  6500   501      - HAVANA       exon    NA    NA   GeneB
#>       gene_name transcript_id exon_number
#> 1 UnsplicedGene          <NA>        <NA>
#> 2          <NA>     GeneA-001        <NA>
#> 3          <NA>     GeneA-001           1
#> 4   SplicedGene          <NA>        <NA>
#> 5          <NA>     GeneB-001        <NA>
#> 6          <NA>     GeneB-001           1

gtf2 <- load_file(tiny_v2_path)
head(gtf2)
#>   seqnames start   end width strand source       type score phase gene_id
#> 1     chr1  1000  5000  4001      + HAVANA       gene    NA    NA   GeneA
#> 2     chr1  1000  5000  4001      + HAVANA transcript    NA    NA   GeneA
#> 3     chr1  1000  5000  4001      + HAVANA       exon    NA    NA   GeneA
#> 4     chr1  6000 12000  6001      - HAVANA       gene    NA    NA   GeneB
#> 5     chr1  6000 12000  6001      - HAVANA transcript    NA    NA   GeneB
#> 6     chr1  6000  6500   501      - HAVANA       exon    NA    NA   GeneB
#>       gene_name transcript_id exon_number
#> 1 UnsplicedGene          <NA>        <NA>
#> 2          <NA>     GeneA-001        <NA>
#> 3          <NA>     GeneA-001           1
#> 4   SplicedGene          <NA>        <NA>
#> 5          <NA>     GeneB-001        <NA>
#> 6          <NA>     GeneB-001           1
```
