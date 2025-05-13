#' @title Compare Annotation Counts Between Two GENCODE Releases
#'
#' @description Compares the number of annotated genomic elements (genes, transcripts, exons, introns) between two specified GENCODE releases. The function accepts input data as either file paths (in GTF/GFF format) or pre-loaded data frames. It computes the absolute difference (delta), the percentage change relative to a chosen baseline, and determines the direction of change (increase or decrease).
#'
#' @usage compare_release(input1, input2, type, gene_type = NULL, baseline = "count2")
#'
#' @param input1 A character string specifying the file path to a GTF/GFF file from the first GENCODE release, or a data frame containing the annotation data.
#' @param input2 A character string specifying the file path to a GTF/GFF file from the second GENCODE release, or a data frame containing the annotation data.
#' @param type A character string indicating the type of genomic element to compare. Valid options are \code{"gene"}, \code{"transcript"}, \code{"exon"}, or \code{"intron"}.
#' @param gene_type An optional character string specifying a particular gene biotype to filter comparisons (e.g., \code{"protein_coding"}, \code{"lncRNA"}). If \code{NULL} (default), all gene types are included.
#' @param baseline A character string defining the baseline for calculating percentage change. Options include:
#' \itemize{
#'   \item \code{"count1"}: Uses the count from the first input (release) as the baseline.
#'   \item \code{"count2"}: Uses the count from the second input (release) as the baseline (default).
#'   \item \code{"average"}: Uses the average of the counts from both inputs as the baseline.
#' }
#'
#' @return A list with the following elements:
#' - `delta`: The absolute difference in the number of annotations.
#' - `percentage`: The percentage change relative to the selected baseline.
#' - `direction`: A string indicating the direction of the change ("increase", "decrease", or "no change").
#'
#' @details
#' This function processes two GENCODE releases to compare annotation counts for a specified genomic element type. Key steps include:
#' \enumerate{
#'   \item Input Handling: If inputs are file paths, they are loaded into data frames using the \code{load_file} function. Data frames are used directly.
#'   \item Element Filtering: If \code{gene_type} is specified, annotations are filtered to include only that gene biotype.
#'   \item Count Calculation: The number of elements (genes, transcripts, etc.) of the specified type is counted in each release.
#'   \item Delta and Percentage: The absolute difference (delta) and percentage change are calculated based on the chosen baseline.
#'   \item Direction Determination: The direction of change is determined by comparing counts between the two releases.
#' }
#' The function provides both numerical results and a formatted console output highlighting key metrics.
#'
#' @examples
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' file_v2 <- system.file("extdata", "gencode.v2.example.gtf.gz", package = "GencoDymo2")
#' # Example 1: Using data frames with the provided example GTF files
#' gtf_v1 <- load_file(file_v1)
#' gtf_v2 <- load_file(file_v2)
#' comparison <- compare_release(gtf_v1, gtf_v2, type = "gene")
#'
#'
#' @seealso \code{\link{load_file}}, \code{\link{get_gtf}}, \code{\link{get_gff3}}.
#' @export

compare_release <- function(input1, input2, type, gene_type = NULL, baseline = "count2") {
  load_if_file <- function(input) {
    if (is.character(input) && file.exists(input)) {
      message("Loading data from file: ", input)
      return(load_file(input)) # Assuming `load_file` function exists
    } else if (is.data.frame(input)) {
      return(input)
    } else {
      stop("Invalid input: Provide either a file path or a data frame.")
    }
  }
  valid_types <- c("gene", "transcript", "exon", "intron")
  if (!type %in% valid_types) {
    stop("Invalid type. Choose from: ", paste(valid_types, collapse = ", "))
  }
  valid_baselines <- c("count1", "count2", "average")
  if (!baseline %in% valid_baselines) {
    stop("Invalid baseline. Choose from: count1, count2, average")
  }
  df1 <- load_if_file(input1)
  df2 <- load_if_file(input2)
  if (!is.null(gene_type)) {
    df1 <- subset(df1, gene_type == gene_type)
    df2 <- subset(df2, gene_type == gene_type)
  }
  count_elements <- function(df, element_type) {
    return(nrow(subset(df, type == element_type)))
  }
  count1 <- count_elements(df1, type)
  count2 <- count_elements(df2, type)
  delta <- abs(count2 - count1)
  if (baseline == "count1") {
    percentage <- round((delta / count1) * 100, digits = 3)
  } else if (baseline == "count2") {
    percentage <- round((delta / count2) * 100, digits = 3)
  } else if (baseline == "average") {
    percentage <- round((delta / ((count1 + count2) / 2)) * 100, digits = 3)
  }
  change_direction <- if (count2 > count1) "increase" else if (count2 < count1) "decrease" else "no change"
  result <- list(delta = delta, percentage = percentage, direction = change_direction)
  cat(paste0("\033[0;32mDelta: \033[0m", delta), sep = "\n")
  cat(paste0("\033[0;32mPercentage: \033[0m", percentage, "%"), sep = "\n")
  cat(paste0("\033[0;32mChange Direction: \033[0m", change_direction), sep = "\n")

  return(invisible(result))
}

#' @title Classify Exons by Their Relative Transcript Position
#'
#' @description Categorizes each exon in a transcript as single, first, inner, or last based on its position. This classification helps in analyzing transcript structure and splicing patterns.
#'
#' @usage classify_exons(input, verbose = TRUE)
#'
#' @param input A character string specifying the path to a GTF/GFF3 file or a data frame containing GTF-formatted data. The file should follow the GENCODE format conventions.
#' @param verbose A logical indicating whether to display progress messages and exon counts. Defaults to \code{TRUE}.
#' @return A data frame identical to the input but with an additional column \code{EXON_CLASSIFICATION}. This column contains one of the following values for each exon:
#' \itemize{
#'   \item \code{"single_exon"}: The transcript contains only one exon.
#'   \item \code{"first_exon"}: The first exon in a multi-exon transcript.
#'   \item \code{"last_exon"}: The last exon in a multi-exon transcript.
#'   \item \code{"inner_exon"}: Exons between the first and last in multi-exon transcripts.
#' }
#'
#' @details The function processes the input GTF data to:
#' \enumerate{
#'   \item Load the GTF file (if a path is provided) or use the provided data frame.
#'   \item Filter entries to include only exons.
#'   \item For each transcript, count the total exons and determine each exon's position.
#'   \item Classify exons based on their position and the total count.
#' }
#' If \code{verbose = TRUE}, the function prints counts of each exon type.
#'
#' @examples
#'
#' # Example 1: Using the provided example GTF files
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' gtf_v1 <- classify_exons(file_v1)
#'
#' # Example 2: Using a pre-loaded data frame
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' gtf_v1 <- load_file(file_v1)
#' classified_data_v1 <- classify_exons(gtf_v1)
#'
#' @import dplyr
#' @export

classify_exons <- function(input, verbose = TRUE) {
  if (is.character(input)) {
    if (verbose) message("Loading GTF file: ", input)
    gtf <- load_file(input)
  } else if (is.data.frame(input)) {
    if (verbose) message("Using input dataframe...")
    gtf <- input
  } else {
    stop("Input must be a file path (character) or a dataframe.")
  }
  exons <- subset(gtf, gtf$type == "exon")
  if (nrow(exons) == 0) {
    stop("Error: No exons found in the input GTF file or dataframe.")
  }
  # Count and classify exons
  if (verbose) message("Classifying exons...")
  exon_classification <- exons %>%
    group_by(transcript_id) %>%
    mutate(
      exon_count = n(),
      EXON_CLASSIFICATION = case_when(
        exon_count == 1 ~ "single_exon",
        exon_number == 1 ~ "first_exon",
        exon_number == exon_count ~ "last_exon",
        TRUE ~ "inner_exon"
      )
    ) %>%
    ungroup() %>%
    select(-exon_count)

  gtf_class <- suppressWarnings(dplyr::left_join(gtf, exon_classification, by = colnames(gtf)))
  single_exons_count <- sum(exon_classification$EXON_CLASSIFICATION == "single_exon")
  first_exons_count <- sum(exon_classification$EXON_CLASSIFICATION == "first_exon")
  last_exons_count <- sum(exon_classification$EXON_CLASSIFICATION == "last_exon")
  inner_exons_count <- sum(exon_classification$EXON_CLASSIFICATION == "inner_exon")
  if (verbose) {
    message(paste0("\033[0;32mSingle exons: \033[0m", single_exons_count))
    message(paste0("\033[0;32mFirst exons: \033[0m", first_exons_count))
    message(paste0("\033[0;32mLast exons: \033[0m", last_exons_count))
    message(paste0("\033[0;32mInner exons: \033[0m", inner_exons_count))
  }
  return(gtf_class)
}

#' @title Eliminate Redundant Genomic Elements
#' @description Removes redundant genomic elements (exons or introns) from a data frame, ensuring each element is uniquely represented. Redundancy is determined by genomic coordinates and gene ID.
#'
#' @usage eliminate_redundant_elements(input, element_type = "exon")
#'
#' @param input A data frame containing genomic element coordinates (exons or introns) with columns \code{seqnames}, \code{start}, \code{end}, and \code{gene_id}.
#' @param element_type The type of genomic element to process. Valid options are \code{"exon"} (default) or \code{"intron"}.
#'
#' @return A data frame with redundant elements removed, retaining only unique entries based on coordinates and gene ID.
#'
#' @details This function uses genomic coordinates (chromosome, start, end) and gene ID to identify and remove duplicate entries. For exons, coordinates are directly compared. For introns, coordinates are derived from \code{intron_start} and \code{intron_end} columns (check \code{extract_introns} function for more details)
#'
#' @examples
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' gtf_v1 <- load_file(file_v1)
#' # Eliminate redundant exons
#' nonredundant_exons <- eliminate_redundant_elements(gtf_v1, element_type = "exon")
#'
#' @importFrom tidyr unite
#' @importFrom data.table rleidv
#' @export

eliminate_redundant_elements <- function(input, element_type = "exon") {
  if (!element_type %in% c("exon", "intron")) {
    stop("Invalid element_type. Please choose either 'exon' or 'intron'.")
  }
  coordinate_columns <- if (element_type == "exon") {
    c("seqnames", "start", "end")
  } else {
    c("seqnames", "intron_start", "intron_end")
  }
  missing_columns <- setdiff(coordinate_columns, colnames(input))
  if (length(missing_columns) > 0) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  uniqueid <- tidyr::unite(input, id, all_of(coordinate_columns), remove = FALSE, sep = "-")
  nonredundant <- uniqueid[!duplicated(data.table::rleidv(uniqueid, cols = c("id", "gene_id"))), ]
  nonredundant$id <- NULL
  return(as.data.frame(nonredundant))
}


#' @title Extract Genomic Elements by Strand
#'
#' @description Filters a data frame or GTF/GFF3 file to extract specific genomic elements (genes, transcripts, exons, or introns) located on a specified DNA strand (+ or -).
#'
#' @usage extract_element_by_strand(input, type, strand, verbose = TRUE)
#'
#' @param input A data frame derived from a GTF/GFF3 file or containing genomic annotations.
#' @param type A character string specifying the type of element to extract. Valid options: \code{"gene"}, \code{"transcript"}, \code{"exon"}, \code{"intron"}.
#' @param strand A character string indicating the DNA strand to filter. Valid options: \code{"+"} (forward) or \code{"-"} (reverse).
#' @param verbose A logical indicating whether to print the count of extracted elements. Defaults to \code{TRUE}.
#'
#' @return A data frame containing only the elements of the specified type located on the chosen strand.
#'
#' @details This function filters the input data based on the \code{type} and \code{strand} columns. It is useful for strand-specific analyses, such as studying antisense transcripts or strand-biased genomic features.
#'
#' @examples
#' # Example 1: Extract genes on the forward strand using gencode.v1
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' gtf_v1 <- load_file(file_v1)
#' forward_genes_v1 <- extract_element_by_strand(gtf_v1, type = "gene", strand = "+")
#'
#' # Example 2: Extract exons on the reverse strand using gencode.v1
#' reverse_exons_v1 <- extract_element_by_strand(gtf_v1, type = "exon", strand = "-")
#'
#' @importFrom dplyr filter
#' @importFrom tools toTitleCase
#' @export

extract_element_by_strand <- function(input, type, strand, verbose = TRUE) {
  if (!is.data.frame(input)) {
    stop("'input' must be a dataframe. Ensure you provide a valid dataframe or GTF data as input.")
  }
  valid_types <- c("gene", "transcript", "exon", "intron")
  if (!type %in% valid_types) {
    stop("Invalid 'type' specified. Choose from: 'gene', 'transcript', 'exon', 'intron'.")
  }
  if (!strand %in% c("+", "-")) {
    stop("Invalid 'strand' specified. Choose from: '+' or '-'.")
  }
  if (strand == "+") {
    filtered_data <- filter(input, type == type & strand == "+")
  } else {
    filtered_data <- filter(input, type == type & strand == "-")
  }
  element_count <- nrow(filtered_data)
  if (verbose) {
    cat(paste0("\033[0;32m", tools::toTitleCase(type), "s on strand '", strand, "': \033[0m", element_count), sep = "\n")
  }

  return(as.data.frame(filtered_data))
}


#' @title Calculate Spliced Transcript Lengths
#'
#' @description Computes the spliced length of transcripts by summing the lengths of their constituent exons. This represents the mature RNA length after intron removal.
#'
#' @usage spliced_trans_length(input)
#'
#' @param input A data frame containing GENCODE annotations in GTF format, including \code{type}, \code{transcript_id}, and \code{width} columns.
#'
#' @return
#' A data frame with two columns: \code{transcript_id} and \code{transcript_length}, where the latter is the sum of exon widths for each transcript.
#'
#' @details This function processes the input data to:
#' \enumerate{
#'   \item Filter entries to include only exons.
#'   \item Group exons by transcript ID.
#'   \item Sum the widths of all exons per transcript.
#' }
#' The result provides the total length of the mature spliced RNA for each transcript.
#'
#' @examples
#'
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' gtf_v1 <- load_file(file_v1)
#' spliced_lengths <- spliced_trans_length(gtf_v1)
#'
#' @importFrom data.table setDT
#' @export

spliced_trans_length <- function(input) {
  if (!is.data.frame(input)) {
    stop("Input must be a dataframe. Ensure you have loaded a GTF file into a dataframe.")
  }
  required_columns <- c("type", "transcript_id", "width")
  missing_columns <- setdiff(required_columns, colnames(input))
  if (length(missing_columns) > 0) {
    stop("Input dataframe is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  exons <- subset(input, type == "exon")
  transcript_lengths <- data.table::setDT(exons)[, .(transcript_length = sum(width)), by = transcript_id]
  result <- as.data.frame(transcript_lengths)
  return(result)
}


#' @title Generate Summary Statistics for Genomic Elements
#'
#' @description Calculates descriptive summary statistics (mean, median, standard deviation, etc.) for the lengths of exons or introns, grouped by classification (e.g., first exon, inner intron).
#'
#' @usage stat_summary(input, type, verbose = TRUE)
#'
#' @param input A data frame containing classified exons or introns. For exons, must include \code{EXON_CLASSIFICATION}. For introns, requires classification columns.
#' @param type A character string specifying the element type. Valid options: \code{"exon"} or \code{"intron"}.
#' @param verbose A logical indicating whether to print progress messages. Defaults to \code{TRUE}.
#'
#' @return A data frame with summary statistics for each element group, including mean, median, standard deviation, standard error, quartiles, and sample size.
#'
#' @details
#' For exons, statistics are grouped by \code{EXON_CLASSIFICATION}. For introns, groups include \code{first_intron}, \code{inner_intron}, and splice site types (e.g., \code{gc_intron}).
#'
#' @examples
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' gtf_v1 <- load_file(file_v1)
#' # Exon statistics
#' exon_stats <- stat_summary(gtf_v1, type = "exon")
#'
#'
#' @import dplyr
#' @importFrom plotrix std.error
#' @export

stat_summary <- function(input, type, verbose = TRUE) {
  if (!type %in% c("exon", "intron")) {
    stop("Invalid type specified. Choose 'exon' or 'intron'.")
  }

  if (type == "exon") {
    if (!"EXON_CLASSIFICATION" %in% colnames(input)) {
      if (verbose) message("Classifying exons...")
      input <- classify_exons(input, verbose = FALSE)
      if (verbose) message("Exon classification completed.")
    } else {
      if (verbose) message("Exons are already classified. Skipping classification.")
    }
    exons <- subset(input, type == "exon")
    data <- dplyr::select(exons, width, group = EXON_CLASSIFICATION)
  } else { # type == "intron"
    if (verbose) message("Classifying introns based on their properties...")
    data_list <- list(
      dplyr::mutate(dplyr::filter(input, intron_number == "intron1"), group = "first_intron"),
      dplyr::mutate(dplyr::filter(input, intron_number != "intron1"), group = "inner_intron")
    )
    data <- dplyr::bind_rows(data_list)

    if (verbose) message("Intron classification completed.")
  }

  if (verbose) message("Calculating summary statistics...")

  # Calculate summary stats and allow NA for empty groups
  stats <- tryCatch(
    {
      data %>%
        dplyr::group_by(group) %>%
        dplyr::summarise(
          mean = mean(width, na.rm = TRUE),
          median = median(width, na.rm = TRUE),
          Q1 = quantile(width, 0.25, na.rm = TRUE),
          Q3 = quantile(width, 0.75, na.rm = TRUE),
          sd = sd(width, na.rm = TRUE),
          std_error = plotrix::std.error(width, na.rm = TRUE),
          N = dplyr::n(),
          .groups = "drop"
        ) %>%
        dplyr::arrange(group)
    },
    error = function(e) {
      warning("Error while calculating summary statistics: ", e$message)
      return(data.frame())
    }
  )

  if (verbose) message("Summary statistics calculation completed.")
  return(as.data.frame(stats))
}


#' @title Calculate GC Content of Genomic Features
#'
#' @description Computes the GC content percentage for exons or introns using genome sequence data from a BSgenome object.
#'
#' @usage calculate_gc_content(input, genome, verbose = TRUE)
#'
#' @param input A data frame containing genomic features (exons/introns) with \code{seqnames}, \code{start}, \code{end}, and \code{strand} columns.
#' @param genome A BSgenome object representing the reference genome (e.g., \code{BSgenome.Hsapiens.UCSC.hg38}).
#' @param verbose A logical indicating whether to print progress messages. Defaults to \code{TRUE}.
#'
#' @return The input data frame with an additional \code{gc_content} column containing GC percentages for each feature.
#'
#' @details This function extracts the DNA sequence for each feature using genomic coordinates, then calculates the proportion of G and C nucleotides. Requires a BSgenome package for the relevant genome.
#'
#' @examples
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' gtf_v1 <- load_file(file_v1)
#' gtf_with_gc <- calculate_gc_content(gtf_v1, genome = BSgenome.Hsapiens.UCSC.hg38, verbose = FALSE)
#'
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings letterFrequency
#' @importFrom GenomicRanges GRanges
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @export

calculate_gc_content <- function(input, genome, verbose = TRUE) {
  if (!methods::is(genome, "BSgenome")) {
    stop("The 'genome' argument must be a BSgenome object.")
  }

  if (!is.data.frame(input)) {
    stop("Input must be a data frame.")
  }

  required_cols <- c("seqnames", "strand", "type")
  if (!all(required_cols %in% colnames(input))) {
    stop("Input data frame must contain at least: ", paste(required_cols, collapse = ", "))
  }

  # Determine if exon or intron coordinates are used
  if (!("start" %in% colnames(input)) && !("intron_start" %in% colnames(input))) {
    stop("Input must contain either 'start' or 'intron_start' column.")
  }
  if (!("end" %in% colnames(input)) && !("intron_end" %in% colnames(input))) {
    stop("Input must contain either 'end' or 'intron_end' column.")
  }

  features <- input[input$type %in% c("exon", "intron"), ]

  if (nrow(features) == 0) {
    if (verbose) message("No exons or introns found in the input data frame. Returning input unchanged.")
    input$gc_content <- NA
    return(input)
  }

  # Add universal start and end for GRanges construction
  features$start_unified <- if ("start" %in% colnames(features)) features$start else features$intron_start
  features$end_unified <- if ("end" %in% colnames(features)) features$end else features$intron_end

  if (verbose) message("Creating GRanges object...")
  gr <- GenomicRanges::GRanges(
    seqnames = features$seqnames,
    ranges = IRanges::IRanges(start = features$start_unified, end = features$end_unified),
    strand = features$strand
  )

  if (verbose) message("Extracting sequences...")
  seqs <- tryCatch(
    {
      BSgenome::getSeq(genome, gr)
    },
    error = function(e) {
      stop("Error extracting sequences: ", e$message)
    }
  )

  if (verbose) message("Calculating GC content...")
  gc_content <- Biostrings::letterFrequency(seqs, letters = "GC", as.prob = TRUE)[, 1]

  input$gc_content <- NA
  input[input$type %in% c("exon", "intron"), "gc_content"] <- gc_content

  if (verbose) message("GC content calculation complete.")

  return(input)
}

#' @title Convert Data Frame to FASTA File
#'
#' @description Converts a data frame containing sequence IDs and sequences into a FASTA-formatted file, optionally compressed as gzip.
#'
#' @usage df_to_fasta(df, id_col, seq_col, output_file = NULL, gzip = TRUE, verbose = TRUE)
#' @param df A data frame with at least two columns: one for sequence IDs and one for sequences.
#' @param id_col A character string specifying the column name containing sequence IDs.
#' @param seq_col A character string specifying the column name containing sequence data.
#' @param output_file A character string specifying the output file path. If NULL, the function will stop with an informative message.
#' @param gzip A logical indicating whether to compress the output as a gzip file. Defaults to \code{TRUE}.
#' @param verbose A logical indicating whether to print progress messages. Defaults to \code{TRUE}.
#'
#' @return No return value. Writes a FASTA file to the specified path.
#'
#' @details This function efficiently writes large sequence datasets to FASTA format, handling compression and progress reporting. It validates input columns and manages memory by processing data in chunks.
#'
#' @examples
#' temp_dir <- tempdir()
#' temp_output <- file.path(temp_dir, "output.fa.gz")
#' seq_data <- data.frame(
#'   transcript_id = c("ENST0001", "ENST0002"),
#'   sequence = c("ATGCTAGCTAG", "GCTAGCTAGCT")
#' )
#' df_to_fasta(seq_data, "transcript_id", "sequence", temp_output)
#'
#' @importFrom progress progress_bar
#' @export

df_to_fasta <- function(df, id_col, seq_col, output_file = NULL, gzip = TRUE, verbose = TRUE) {
  if (is.null(output_file)) {
    stop("The 'output_file' argument cannot be NULL. Please specify a valid output path.")
  }
  if (!id_col %in% colnames(df)) {
    stop(paste("The column", id_col, "is not present in the dataframe."))
  }
  if (!seq_col %in% colnames(df)) {
    stop(paste("The column", seq_col, "is not present in the dataframe."))
  }
  if (gzip && !grepl("\\.gz$", output_file)) {
    output_file <- paste0(output_file, ".gz")
  } else if (!gzip && grepl("\\.gz$", output_file)) {
    output_file <- sub("\\.gz$", "", output_file)
  }
  fasta_entries <- paste0(">", df[[id_col]], "\n", df[[seq_col]])
  con <- if (gzip) gzfile(output_file, "w") else file(output_file, "w")
  if (verbose) {
    cat("Writing FASTA file...\n")
    pb <- progress::progress_bar$new(
      format = "[:bar] :current/:total (:percent) in :elapsed | ETA: :eta",
      total = length(fasta_entries),
      clear = FALSE,
      width = 60
    )
  }
  chunk_size <- 10000
  for (i in seq(1, length(fasta_entries), by = chunk_size)) {
    chunk <- fasta_entries[i:min(i + chunk_size - 1, length(fasta_entries))]
    cat(chunk, file = con, sep = "\n")
    if (verbose) {
      pb$tick(length(chunk))
    }
  }
  close(con)
  if (verbose) {
    cat("FASTA file successfully saved to", output_file, "\n")
  }
}



#' @title Extract Coding Sequences (CDS) from GTF Annotations
#'
#' @description Extracts CDS regions from a GTF annotation file or data frame using genomic coordinates and retrieves corresponding DNA sequences from a BSgenome reference.
#'
#' @usage extract_cds_sequences(input, genome, save_fasta, output_file, verbose)
#'
#' @param input A character string (GTF file path) or data frame containing CDS annotations.
#' @param genome A BSgenome object for the relevant genome. Defaults to human (hg38).
#' @param save_fasta A logical indicating whether to save sequences to a FASTA file. Defaults to \code{FALSE}.
#' @param output_file A character string specifying the FASTA output path. If \code{NULL}, uses "CDS.fa".
#' @param verbose A logical indicating whether to print progress messages. Defaults to \code{TRUE}.
#'
#' @return A data frame containing CDS annotations with corresponding sequences. If \code{save_fasta = TRUE}, also writes a FASTA file.
#'
#' @details This function processes CDS entries from the input GTF, extracts their sequences from the reference genome, and optionally saves them in FASTA format. Useful for downstream analyses like protein translation.
#'
#' @examples
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' gtf_v1 <- load_file(file_v1)
#' # Human CDS extraction
#' suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
#' suppressPackageStartupMessages(library(GenomicRanges))
#' gtf_granges <- GRanges(gtf_v1)
#' cds_seqs <- extract_cds_sequences(gtf_granges, BSgenome.Hsapiens.UCSC.hg38, save_fasta = FALSE)
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings DNAStringSet
#' @importFrom BSgenome getSeq
#' @importFrom methods is
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @export

extract_cds_sequences <- function(input, genome = BSgenome.Hsapiens.UCSC.hg38, save_fasta = FALSE, output_file = NULL, verbose = TRUE) {
  if (!methods::is(genome, "BSgenome")) {
    stop("The 'genome' argument must be a BSgenome object.")
  }
  # Handle input
  if (is.character(input)) {
    if (!file.exists(input)) {
      stop("GTF file not found: ", input)
    }
    if (verbose) message("Loading GTF file...")
    gtf <- load_file(input) # returns GRanges
  } else if (methods::is(input, "GRanges")) {
    if (verbose) message("Using provided GRanges object...")
    gtf <- input
  } else {
    stop("Input must be a file path or a GRanges object.")
  }
  gtf <- GenomicRanges::GRanges(gtf)

  # Extract CDS features
  cds <- subset(gtf, gtf$type == "CDS")

  if (length(cds) == 0) {
    warning("No CDS features found in the GTF data.")
    return(Biostrings::DNAStringSet())
  }

  if (verbose) message("Extracting CDS sequences from genome...")

  cds_seqs <- tryCatch(
    {
      BSgenome::getSeq(genome, cds)
    },
    error = function(e) {
      stop("Error extracting sequences from genome: ", e$message)
    }
  )

  # Prepare output table
  cds_df <- as.data.frame(cds)
  cds_df$CDS_seq <- as.character(cds_seqs)

  # Optionally write FASTA
  if (save_fasta) {
    if (verbose) cat("Compiling sequences into a FASTA file...\n")
    if (is.null(output_file)) {
      output_file <- "CDS.fa"
    }
    cds_df$ID <- paste0(cds_df$transcript_id, ";", cds_df$exon_id)
    df_to_fasta(cds_df, id_col = "ID", seq_col = "CDS_seq", output_file = output_file, verbose = TRUE)
    cds_df$ID <- NULL
  }
  return(cds_df)
}
