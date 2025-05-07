#' @name tiny_example_gtf_files
#' @title Tiny example GTF files
#'
#' @description
#' These are minimal GTF files used for examples and testing within the package.
#'
#' - **gencode.v1.example.gtf.gz** contains two genes:
#'   - GeneA: a single-exon, unspliced gene.
#'   - GeneB: a spliced gene with two transcripts and multiple exons.
#'
#' - **gencode.v2.example.gtf.gz** contains the same two genes as in `gencode.v1.example.gtf.gz` plus:
#'   - GeneC: a new spliced gene with multiple transcripts and many exons.
#'
#' These files are stored in the `inst/extdata/` directory and can be accessed using
#' \code{system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")}
#' or
#' \code{system.file("extdata", "gencode.v2.example.gtf.gz", package = "GencoDymo2")}.
#'
#' @format Two external GTF files.
#' @seealso \code{\link{load_file}}
#' @source Generated manually for testing purposes.
#' @examples
#' tiny_v1_path <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' tiny_v2_path <- system.file("extdata", "gencode.v2.example.gtf.gz", package = "GencoDymo2")
#'
#' gtf1 <- load_file(tiny_v1_path)
#' head(gtf1)
#'
#' gtf2 <- load_file(tiny_v2_path)
#' head(gtf2)
NULL
