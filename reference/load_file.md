# Load a GTF or GFF3 file from GENCODE as a data frame.

This function imports a GTF or GFF3 file (commonly from the GENCODE
website) and converts it into a data frame. The function provides
flexibility for users to work with genomic feature files easily in the R
environment.

## Usage

``` r
load_file(filename)
```

## Arguments

- filename:

  A character string representing the path to the GTF or GFF3 file
  (e.g., "gencode.vM36.annotation.gtf.gz"). The file could be in GTF or
  GFF3 format and must be downloaded from a reliable source like
  GENCODE.

## Value

A data frame containing the parsed content of the GTF or GFF3 file. The
data frame includes standard columns such as 'seqnames', 'start', 'end',
'strand', 'feature', and 'gene_id', among others.

## Details

The function uses the `rtracklayer` package to import the GTF or GFF3
file and returns it as a data frame. The user should ensure that the
input file is properly formatted and accessible from the specified path.
Files larger than a few hundred MBs may take longer to load and process.

## Examples

``` r
# Load example GTF files from the package
file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
gtf_v1 <- load_file(file_v1)
head(gtf_v1)
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
