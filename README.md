
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GencoDymo2

<!-- badges: start -->

[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/monahton)
![GitHub
issues](https://img.shields.io/github/issues/monahton/GencoDymo2)
![GitHub last
commit](https://img.shields.io/github/last-commit/monahton/GencoDymo2)
![Platform](https://img.shields.io/badge/platform-all-green)

[![Last-changedate](https://img.shields.io/badge/last%20change-2021--11--25-yellowgreen.svg)](/commits/master)
<!-- badges: end -->

------------------------------------------------------------------------

# Data Extraction and Manipulation from the GENCODE Database

## Introduction

The **GENCODE** project \[1\], part of the ENCODE project \[2\], offers
accurate annotations of the human and mouse genomes derived from
literature, primary data repositories, computational predictions, manual
annotation, and experimental validation of genes and transcripts,
including noncoding transcripts that are continually discovered and
annotated. GENCODE gene annotations are regularly updated as the
Ensembl/GENCODE gene sets, accessible via the official website
(<https://www.gencodegenes.org>) and are released approximately four
times a year for mouse and twice a year for human \[1\].​

**GencoDymo2** is an R package designed to interrogate different
releases of GENCODE annotations from the human and mouse genomes. It
provides a streamlined interface for accessing and processing GENCODE
annotations, providing easily accessible data regarding annotation
statistics, release comparisons, and information on introns and splice
sites. Moreover, GencoDymo2 can produce FASTA files of donor and
acceptor splice site motifs that can be directly uploaded to the
MaxEntScan tool \[3\] for the calculation of splice site scores.

## Installation

You can install the development version of GencoDymo2 from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("monahton/GencoDymo2")
```

After installation, load the package:

``` r
library(GencoDymo2)
```

## Usage

For details on how to use this package, please see the vignette.  
Additional details on how to install and use GencoDymo2 are available on
the following link:

In case of problems installing or using the package, open an issue or
send an email to <aboualezz.monah@hsr.it>.
