# GencoDymo2: Comprehensive Analysis of GENCODE Annotations and Splice Site Motifs

------------------------------------------------------------------------

  

## Introduction

------------------------------------------------------------------------

Gene annotations are provided by several databases that integrate data
from systematic explorations of various experimental and computational
resources. The **GENCODE** project \[1\], part of the ENCODE project
\[2\], offers accurate annotations of the human and mouse genomes
derived from literature, primary data repositories, computational
predictions, manual annotation, and experimental validation of genes and
transcripts, including noncoding transcripts that are continually
discovered and annotated. GENCODE is one of the most comprehensive and
standardized databases for gene annotations, widely used by the
scientific community.​

The gene sets provided by GENCODE are comprehensive and include
protein-coding and non-coding loci, encompassing alternatively spliced
isoforms and pseudogenes. GENCODE gene annotations are regularly updated
and released as the Ensembl/GENCODE gene sets, accessible via the
official website <https://www.gencodegenes.org>. As of October 2024, the
latest human release is GENCODE 47, and the latest mouse release is
GENCODE M36 \[1\]. GENCODE gene sets are released approximately four
times a year for mouse and twice a year for human \[1\].​

**GencoDymo2** is an R package designed to interrogate different
releases of GENCODE annotations from the human and mouse genomes. It
provides a streamlined interface for accessing and processing GENCODE
annotations, providing easily accessible data regarding annotation
statistics, release comparisons, and information on introns and splice
sites. Moreover, GencoDymo2 can produce FASTA files of donor and
acceptor splice site motifs that can be directly uploaded to the
MaxEntScan tool \[3\] for the calculation of splice site scores.

  

GencoDymo2 can be used to import and process any gtf/gff3 formatted
file. It automatizes and speeds up the following data manipulation
steps:

- Dynamic detection of latest GENCODE releases from the official website
- Unified interface for GTF/GFF3 file retrieval for the human and mouse
  genomes
- Compare annotations from different GENCODE releases
- Perform summary statistics on gene annotations
- Extract introns information
- Assign splice sites consensus sequence
- Compile MaxEntScan splice sites motifs

------------------------------------------------------------------------

  

## Installation

------------------------------------------------------------------------

GencoDymo2 runs in the R statistical computing environment. You will
need R version 4.1.0 or higher. The most recent version of GencoDymo2
can be installed from GitHub. We recommend using the pak package for
fast and reliable installation.

``` r
# Install pak if not already installed
if (!require("pak")) install.packages("pak")
# Install from GitHub
pak::pkg_install("github::monahton/GencoDymo2")

# Load the package
library(GencoDymo2)
```

------------------------------------------------------------------------

  

## Usage

------------------------------------------------------------------------

### Latest Releases

GENCODE releases annotations approximately every three months,
coinciding with Ensembl releases. The latest release from GENCODE for
human and mouse genomes can be obtained using the
[`get_latest_release()`](https://monahton.github.io/GencoDymo2/reference/get_latest_release.md)
function. This function queries the official GENCODE website to
determine the most recent release for a given species (human or mouse).

``` r
# Fetch the most recent human and mouse GENCODE release identifiers
human_release <- get_latest_release("human", verbose = T)
mouse_release <- get_latest_release("mouse", verbose = T)
```

    ## Latest human GENCODE release: release_47

    ## Latest human GENCODE release: release_M36

  

> ***NOTE:*** Behind the scenes, get_latest_release() parses the Ensembl
> public FTP site, extracts the release version (e.g., “release_47”),
> and returns it for programmatic use.

  

### Download GTF/GFF3 Files

GencoDymo2 allows the download of `gtf` and `gff3` files of the latest
release or any release specified by the user. These files can be
obtained for human and mouse genome and can be saved in a path specified
by the user. With the release identifiers in hand, we can download these
annotation files:

``` r
# Download latest human long noncoding RNAs GTF
lnc_47_gtf <- get_gtf(
  species = "human",
  release_version = human_release,
  annotation_type = "long_noncoding_RNAs.gtf.gz",
  dest_folder = tempdir()
)

# Download previous human release (release_46) for comparison
lnc_46_gtf <- get_gtf(
  species = "human",
  release_version = "release_46",
  annotation_type = "long_noncoding_RNAs.gtf.gz",
  dest_folder = tempdir()
)

# Download latest mouse primary assembly annotations (GFF3)
mouse_36_gff3 <- get_gff3(
  species = "mouse",
  release_version = mouse_release,
  annotation_type = "primary_assembly.annotation.gff3.gz",
  dest_folder = tempdir()
)
```

#### Supported Annotation Types

The valid annotation types for both the gtf and gff3 format can be one
of the following:

    ## Valid Annotation Types:

    ## [1] "annotation"                           
    ## [2] "basic.annotation"                     
    ## [3] "chr_patch_hapl_scaff.annotation"      
    ## [4] "chr_patch_hapl_scaff.basic.annotation"
    ## [5] "long_noncoding_RNAs"                  
    ## [6] "primary_assembly.annotation"          
    ## [7] "primary_assembly.basic.annotation"    
    ## [8] "tRNAs"                                
    ## [9] "polyAs"

  

### Load and Inspect

A `gtf` or `gff3` file can be imported to R as a dataframe using the
[`load_file()`](https://monahton.github.io/GencoDymo2/reference/load_file.md)
function. Once downloaded, load the files into R as data frames using
load_file(). This function automatically detects GTF vs GFF3 and parses
the attributes into columns.

``` r
# Loading using the stored paths from previous steps
lnc_47_df <- load_file(lnc_47_gtf)
head(lnc_47_df)

# Alternatively, specify the file path directly
lnc_46_df <- load_file(file.path(tempdir(), "gencode.v46.long_noncoding_RNAs.gtf.gz"))
head(lnc_46_df)

# Load mouse GFF3
mouse_pri_36 <- load_file(file.path(tempdir(),"gencode.vM36.primary_assembly.annotation.gff3.gz"))
head(mouse_pri_36)
```

  

### Compare Releases

As GENCODE regularly updates its annotations, different releases from
GENCODE contain differing number of annotations. Some genes are added in
new releases while others are classified into different gene classes. In
the majority of cases, the Ensembl gene ID of genes remains constant
while the gene name might vary.  
To quantify changes between two releases, use
[`compare_release()`](https://monahton.github.io/GencoDymo2/reference/compare_release.md).
Specify the two annotation data frames and the type of feature to
compare (“gene”, “transcript”, “exon”, etc.). Additional filters like
the gene_biotype and the baseline can allow to focus on subsets (e.g.,
only TEC or protein_coding genes) and different reference counts for the
calculation of the difference percentage.

``` r
# Compare gene counts between release 47 and 46
gene_comparison <- compare_release(lnc_47_df, lnc_46_df, type = "gene")

# Compare exon counts
exon_comparison <- compare_release(lnc_47_df, lnc_46_df, type = "exon")

# Compare a specific gene biotype (e.g., TEC) using a custom baseline
comparison <- compare_release(
  lnc_47_df,
  lnc_46_df,
  type = "gene",
  gene_biotype = "TEC",
  baseline = "count1"
)
```

The function produces as an output a ***Delta*** value, that corresponds
to the difference in number between the specified elements of the two
inputs, in addition to its ***Percentage***.

  

> ***NOTE:*** The output is a data frame summarizing absolute and
> percentage differences in feature counts, enabling rapid assessment of
> annotation growth or refinement.

  

### Extract Introns Coordinates

The gtf files provided by GENCODE list genes, transcripts, and exons
while intronic regions must be inferred. Given the loaded data frame
from the gtf files of GENCODE, the
[`extract_introns()`](https://monahton.github.io/GencoDymo2/reference/extract_introns.md)
function computes intron coordinates for each transcript by sorting exon
ranges and identifying gaps. It returns a dataframe containing introns
coordinates and position within each transcript.

``` r
# Human lncRNA introns for release 47
introns_lnc_47 <- extract_introns(lnc_47_df, verbose = T)

# Mouse introns (filtering to primary chromosomes first)
mouse_pri_36 <- mouse_pri_36[grepl("^chr", mouse_pri_36$seqnames), ]
mouse_introns_pri_36 <- extract_introns(mouse_pri_36, verbose = T)
```

  

### Extract Donor and Acceptor Splice Sites

Splice sites are defined by dinucleotide motifs at intron boundaries. At
the 5’ end, the donor site includes an almost invariant sequence GT (GU
in the RNA). At the 3’ end, the intron terminates with the acceptor
site, an almost invariant AG sequence.

  

![Fig1. Donor and acceptor splice sites consensus. Image taken from
(https://commons.wikimedia.org/wiki/File:Consensus_intron.jpg).](Consensus_intron.jpg)

Fig1. Donor and acceptor splice sites consensus. Image taken from
(<https://commons.wikimedia.org/wiki/File:Consensus_intron.jpg>).

  

The majority of introns belong to the U2-type spliceosome and are
flanked by GT– AG splice site dinucleotides. The most frequent exception
to this rule are the U2-type GC–AG splice sites, comprising ∼0.8% of
human splice sites and about ∼0.09% of the human splice sites belong to
the U12-type which are processed by the minor spliceosome and are
flanked by AT–AC dinucleotides \[4\].  
The
[`assign_splice_sites()`](https://monahton.github.io/GencoDymo2/reference/assign_splice_sites.md)
function retrieves these flanking sequences from a BSgenome object. It
adds two new columns to the data frame containing the introns
coordinates, assigning both the 5’ and 3’ splice sites. Sequences should
be retrieved from a loaded BSgenome object such as
`BSgenome.Hsapiens.UCSC.hg38` for human and
`BSgenome.Mmusculus.UCSC.mm39` for mouse.

``` r
# Human
library(BSgenome.Hsapiens.UCSC.hg38)
lnc_47_ss <- assign_splice_sites(
  introns_lnc_47,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  verbose = T
)

# Mouse
library(BSgenome.Mmusculus.UCSC.mm39)
mouse_pri_36_ss <- assign_splice_sites(
  mouse_introns_pri_36,
  genome = BSgenome.Mmusculus.UCSC.mm39,
  verbose = T
)
```

  

#### Identify Cryptic splice sites

Cryptic splice sites are those that do not match the canonical donor
(GT) or acceptor motifs (AG). The
[`find_cryptic_splice_sites()`](https://monahton.github.io/GencoDymo2/reference/find_cryptic_splice_sites.md)
function identifies potential cryptic splice sites by comparing sequence
motifs in introns to canonical splice site motifs (donor and acceptor)
that can be also specified by the user. It compares the identified
splice sites with the provided canonical motifs and flags the sites that
differ from the canonical patterns, making it useful for studying
aberrant splicing events.

``` r
# Identify cryptic (non-canonical) splice sites
cryptic_ss <- find_cryptic_splice_sites(
  lnc_47_ss,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  canonical_donor = "GT",
  canonical_acceptor = "AG",
  verbose = TRUE
)
```

  

#### Extract splice sites motifs for MaxEntScan webtool

**MaxEntScan** is a webtool based on the approach for modeling the
sequences of short sequence motifs involved in RNA splicing which
simultaneously accounts for non-adjacent as well as adjacent
dependencies between positions. This method is based on the *‘Maximum
Entropy Principle’* and generalizes most previous probabilistic models
of sequence motifs such as weight matrix models and inhomogeneous Markov
models \[3\].  
The **MaxEntScan::score5ss** assign scores for donor splice sites
according to four models by using **9-mer SEQUENCES** motif as input in
a `.fa file`. The 9-mers motif contain 3 bases in the exon and 6 bases
in the intron.  
The **MaxEntScan::score3ss** assign scores for donor splice sites
according to three models. It takes as input a `.fa file` of **23-mer
SEQUENCES** in which 20 bases are in the intron and 3 bases in the
exon.  
The
[`extract_ss_motif()`](https://monahton.github.io/GencoDymo2/reference/extract_ss_motif.md)
function is used along with BSgenome object of the studied species to
retrieve *MaxEntScan::score5ss* and *MaxEntScan::score3ss* motif
sequences.  
The user can choose to generate a **fasta file** in the working
directory or any specified path as an output, that contains either
9-mers or 23-mers, respectively for each 5’ and 3’ splice-sites. It
generates also a **dataframe** for the coordinates and IDs of the
corresponding motifs.  
The generated fasta files can then be directly utilized by *MaxEntScan*
tools.  
The first argument `input` is a dataframe containing the intron
coordinates. The second argument `genome` is a BSgenome object (the
human genome sequence hg38 is used by default). The third argument
`type` specifies whether to extract motifs at the 5ss or the 3ss. Users
can optionally set the argument `save_fasta` to TRUE to and the
`output_file` argument to save a fasta file.

``` r
# Donor motifs (5'ss)
motifs_donor <- extract_ss_motif(
  input = lnc_47_ss,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  type = "5ss",
  verbose = T,
  save_fasta = T,
  output_file = file.path(tempdir(), "lnc_47_5ss_motifs.fa")
)

# Acceptor motifs (3'ss)
motifs_acc <- extract_ss_motif(
  input = lnc_47_ss,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  type = "3ss",
  verbose = T,
  save_fasta = T,
  output_file = file.path(tempdir(), "lnc_47_3ss_motifs.fa")
)
```

------------------------------------------------------------------------

  

## Additional Utilities

------------------------------------------------------------------------

#### Single-Exon Features

Some genes are unspliced: composed by one single exon and have no
introns.  
The
[`extract_single_exon()`](https://monahton.github.io/GencoDymo2/reference/extract_single_exon.md)
function can extract all the ***single-exon*** genes or transcripts and
returns in a dataframe with their ensembl gene IDs or transcript IDs.

``` r
## identify single exon genes and transcripts
single_exon_genes <- extract_single_exon(lnc_47_df, level = "gene")
single_exon_trans <- extract_single_exon(lnc_47_df, level = "transcript")
```

  

#### Exon Classification

The
[`classify_exons()`](https://monahton.github.io/GencoDymo2/reference/classify_exons.md)
adds a new column `(EXON_CLASSIFICATION)` to the input dataframe,
describing the position of each exon. It labels each exon as
*(first_exons)*, *(inner_exons)*, *(last_exons)* and *(single_exons)*
facilitating position-based analyses. This classification is informative
considering that first and last exons often have peculiar regulative
roles.

``` r
# Assign the ordinal position of exons
lnc_47_class_exons <- classify_exons(lnc_47_df, verbose = TRUE)
```

  

#### Spliced Transcript Length

The
[`spliced_trans_length()`](https://monahton.github.io/GencoDymo2/reference/spliced_trans_length.md)
function computes the mature (post-splicing) transcript length by
summing exon widths per transcript.

``` r
# Length of spliced transcript  
lnc_47_spliced_length <- spliced_trans_length(lnc_47_df)
head(lnc_47_spliced_length)
```

  

#### Summary Statistics

GencoDymo2 provides the
[`stat_summary()`](https://monahton.github.io/GencoDymo2/reference/stat_summary.md)
function that gives a brief summary of data annotations in a particular
GENCODE release regarding genes, transcripts, exons and introns.

``` r
# Exon length statistics
lnc_47_exon_stats <- stat_summary(lnc_47_class_exons, type = "exon")

# Intron length statistics
lnc_47_intron_stats <- stat_summary(introns_lnc_47, type = "intron")
```

  

#### GC Content Calculation

The
[`calculate_gc_content()`](https://monahton.github.io/GencoDymo2/reference/calculate_gc_content.md)
function computes GC percentage for any feature set using a BSgenome
reference.

``` r
# Human
lnc_47_gc <- calculate_gc_content(
  lnc_47_df,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  verbose = TRUE
)
# Mouse
mouse_pri_36_gc <- calculate_gc_content(
  mouse_pri_36,
  genome = BSgenome.Mmusculus.UCSC.mm39,
  verbose = TRUE
)
```

  

#### Extract CDS Sequences

For coding annotations,
[`extract_cds_sequences()`](https://monahton.github.io/GencoDymo2/reference/extract_cds_sequences.md)
writes CDS FASTA files from GRanges objects.

``` r
# Convert to GRanges and extract
library(GenomicRanges)
mouse_pri_36_granges <- GRanges(mouse_pri_36)
mouse_cds_seqs <- extract_cds_sequences(
  mouse_pri_36_granges,
  BSgenome.Mmusculus.UCSC.mm39,
  save_fasta = TRUE,
  output_file = file.path(tempdir(), "mouse_pri_36_CDS.fa.gz")
  verbose = TRUE
)
```

------------------------------------------------------------------------

  

## Session Information

------------------------------------------------------------------------

The following provides the session information used when compiling this
document.

    ## ─ Session info ───────────────────────────────────────────────────────────────
    ##  setting  value
    ##  version  R version 4.5.2 (2025-10-31)
    ##  os       Ubuntu 24.04.3 LTS
    ##  system   x86_64, linux-gnu
    ##  ui       X11
    ##  language en-US
    ##  collate  C.UTF-8
    ##  ctype    C.UTF-8
    ##  tz       UTC
    ##  date     2026-01-09
    ##  pandoc   3.1.11 @ /opt/hostedtoolcache/pandoc/3.1.11/x64/ (via rmarkdown)
    ##  quarto   NA
    ## 
    ## ─ Packages ───────────────────────────────────────────────────────────────────
    ##  package     * version date (UTC) lib source
    ##  bslib         0.9.0   2025-01-30 [1] RSPM
    ##  cachem        1.1.0   2024-05-16 [1] RSPM
    ##  cli           3.6.5   2025-04-23 [1] RSPM
    ##  desc          1.4.3   2023-12-10 [1] RSPM
    ##  devtools      2.4.6   2025-10-03 [1] RSPM
    ##  digest        0.6.39  2025-11-19 [1] RSPM
    ##  ellipsis      0.3.2   2021-04-29 [1] RSPM
    ##  evaluate      1.0.5   2025-08-27 [1] RSPM
    ##  fastmap       1.2.0   2024-05-15 [1] RSPM
    ##  fs            1.6.6   2025-04-12 [1] RSPM
    ##  glue          1.8.0   2024-09-30 [1] RSPM
    ##  htmltools     0.5.9   2025-12-04 [1] RSPM
    ##  htmlwidgets   1.6.4   2023-12-06 [1] RSPM
    ##  jquerylib     0.1.4   2021-04-26 [1] RSPM
    ##  jsonlite      2.0.0   2025-03-27 [1] RSPM
    ##  knitr         1.51    2025-12-20 [1] RSPM
    ##  lifecycle     1.0.5   2026-01-08 [1] RSPM
    ##  magrittr      2.0.4   2025-09-12 [1] RSPM
    ##  memoise       2.0.1   2021-11-26 [1] RSPM
    ##  otel          0.2.0   2025-08-29 [1] RSPM
    ##  pkgbuild      1.4.8   2025-05-26 [1] RSPM
    ##  pkgdown       2.2.0   2025-11-06 [1] RSPM
    ##  pkgload       1.4.1   2025-09-23 [1] RSPM
    ##  purrr         1.2.0   2025-11-04 [1] RSPM
    ##  R6            2.6.1   2025-02-15 [1] RSPM
    ##  ragg          1.5.0   2025-09-02 [1] RSPM
    ##  remotes       2.5.0   2024-03-17 [1] RSPM
    ##  rlang         1.1.6   2025-04-11 [1] RSPM
    ##  rmarkdown     2.30    2025-09-28 [1] RSPM
    ##  sass          0.4.10  2025-04-11 [1] RSPM
    ##  sessioninfo   1.2.3   2025-02-05 [1] RSPM
    ##  systemfonts   1.3.1   2025-10-01 [1] RSPM
    ##  textshaping   1.0.4   2025-10-10 [1] RSPM
    ##  usethis       3.2.1   2025-09-06 [1] RSPM
    ##  vctrs         0.6.5   2023-12-01 [1] RSPM
    ##  xfun          0.55    2025-12-16 [1] RSPM
    ##  yaml          2.3.12  2025-12-10 [1] RSPM
    ## 
    ##  [1] /home/runner/work/_temp/Library
    ##  [2] /opt/R/4.5.2/lib/R/site-library
    ##  [3] /opt/R/4.5.2/lib/R/library
    ## 
    ## ──────────────────────────────────────────────────────────────────────────────

------------------------------------------------------------------------

  

## References

------------------------------------------------------------------------

1.  Frankish, A., Diekhans, M., Ferreira, A.M., et al. (2023). *GENCODE:
    reference annotation for the human and mouse genomes in 2023*.
    Nucleic Acids Research, 51(D1), D942–D949.  
2.  ENCODE Project Consortium. (2012). *An integrated encyclopedia of
    DNA elements in the human genome*. Nature, 489(7414), 57–74.  
3.  Yeo, G., & Burge, C.B. (2004). *Maximum entropy modeling of short
    sequence motifs with applications to RNA splicing signals*. Journal
    of Computational Biology, 11(2-3), 377–394.  
4.  Parada GE, Munita R, Cerda CA and Gysling K. *“A comprehensive
    survey of non-canonical splice sites in the human transcriptome”*.
    Nucleic Acids Res. 2014. 42: 10564-10578.  
5.  Abou Alezz, M., Celli, L., Belotti, G., et al. (2020). *GC-AG
    Introns Features in Long Non-coding and Protein-Coding Genes Suggest
    Their Role in Gene Expression Regulation*. Front. Genet. 2020
    11:488.
