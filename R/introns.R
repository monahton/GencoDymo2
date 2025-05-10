#' @title Extract Intron Coordinates from GENCODE Annotations
#'
#' @description Processes a GTF file or data frame to extract intron coordinates, including their genomic positions, transcript associations, and metadata. The function handles both positive and negative strands, ensuring correct orientation of intron boundaries. It is designed to work with GENCODE-formatted annotations.
#'
#' @usage extract_introns(input, verbose = TRUE)
#'
#' @param input A character string specifying the file path to a GTF/GFF3 file, or a data frame containing GTF data. The input must include columns: \code{seqnames}, \code{start}, \code{end}, \code{strand}, \code{type}, and \code{transcript_id}.
#' @param verbose Logical. If \code{TRUE} (default), progress messages are printed during execution.
#'
#' @return A data frame with the following columns:
#' \itemize{
#'   \item \code{seqnames}: Chromosome or scaffold name.
#'   \item \code{intron_start}, \code{intron_end}: Genomic start/end positions of the intron.
#'   \item \code{width}: Length of the intron (intron_end - intron_start + 1).
#'   \item \code{gene_id}, \code{transcript_id}, \code{intron_number}: Gene/transcript identifiers and intron position within the transcript.
#'   \item Additional metadata columns from the original GTF (e.g., \code{gene_name}).
#' }
#'
#' @details
#' \enumerate{
#'   \item Input Handling: Loads the GTF file if a path is provided; uses the data frame directly if supplied.
#'   \item Exon Filtering: Extracts exon records and classifies them (e.g., first/last exon) if needed.
#'   \item Multi-Exon Transcripts: Removes single-exon transcripts to focus on spliced transcripts.
#'   \item Intron Calculation: Determines intron coordinates by identifying gaps between consecutive exons, adjusting for strand orientation.
#'   \item Output: Returns a data frame with intron coordinates and metadata, sorted by gene and position.
#' }
#'
#' @examples
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#'
#' # From a pre-loaded data frame
#' gtf_v1 <- load_file(file_v1)
#' introns <- extract_introns(gtf_v1)
#'
#' @seealso \code{\link{load_file}}, \code{\link{classify_exons}}
#' @import dplyr
#' @export

extract_introns <- function(input, verbose = TRUE) {
  if (is.character(input)) {
    if (!file.exists(input)) {
      stop("File not found: ", input)
    }
    if (verbose) message("Loading GTF file...")
    input <- load_file(input)
  }
  exons <- subset(input, input$type == "exon")

  # Checking if EXON_CLASSIFICATION column exists; if not, classify exons
  if (!"EXON_CLASSIFICATION" %in% colnames(exons)) {
    if (verbose) message("Preparing input...")
    classified_exons <- classify_exons(exons, verbose = FALSE)
  } else {
    if (verbose) message("EXON_CLASSIFICATION already present.")
    classified_exons <- exons
  }

  # Removing single-exon transcripts
  if (verbose) message("Removing single-exon transcripts...")
  transcript_info <- subset(classified_exons, select = c("transcript_id", "exon_number"))
  transcript_counts <- as.data.frame(table(transcript_info$transcript_id))
  single_exon_transcripts <- subset(transcript_counts, Freq == "1")
  if (verbose) message("Single-exon transcripts: ", nrow(single_exon_transcripts))
  multi_exon_trans <- suppressWarnings(
    dplyr::anti_join(classified_exons, single_exon_transcripts, by = c("transcript_id" = "Var1"))
  )
  colnames(transcript_counts) <- c("transcript_id", "exon_count")

  # Extracting intron coordinates
  if (verbose) message("Extracting intron coordinates...")
  exon_count_values <- sort(unique(as.integer(transcript_counts$exon_count)))
  max_exon_count <- max(exon_count_values)
  exon_numbers <- sort(unique(as.integer(multi_exon_trans$exon_number)))
  intersecting_exons <- dplyr::setdiff(
    union(exon_count_values, exon_numbers),
    intersect(exon_count_values, exon_numbers)
  )

  ## Last introns
  last_introns <- NULL
  for (exon_index in exon_count_values) {
    if (exon_index == max_exon_count) {
      next
    }
    current_exons <- subset(
      multi_exon_trans,
      multi_exon_trans$EXON_CLASSIFICATION != "last_exon" & multi_exon_trans$exon_number == exon_index
    )
    current_exons$intron_start <- ifelse(
      current_exons$strand == "+",
      current_exons$end + 1,
      current_exons$start - 1
    )
    next_exons <- subset(multi_exon_trans, multi_exon_trans$exon_number == exon_index + 1)
    next_exons$intron_end <- ifelse(
      next_exons$strand == "+",
      next_exons$start - 1,
      next_exons$end + 1
    )
    intron_ends <- subset(next_exons, select = "intron_end")
    intron_data <- cbind(current_exons, intron_ends)
    intron_data$intron_number <- paste("intron", exon_index, sep = "")
    last_introns <- rbind(last_introns, intron_data)
  }

  ## Inner introns
  inner_introns <- NULL
  for (exon_index in intersecting_exons) {
    current_exons <- subset(multi_exon_trans, multi_exon_trans$exon_number == exon_index)
    current_exons$intron_start <- ifelse(
      current_exons$strand == "+",
      current_exons$end + 1,
      current_exons$start - 1
    )
    next_exons <- subset(multi_exon_trans, multi_exon_trans$exon_number == exon_index + 1)
    next_exons$intron_end <- ifelse(
      next_exons$strand == "+",
      next_exons$start - 1,
      next_exons$end + 1
    )
    intron_ends <- subset(next_exons, select = "intron_end")
    intron_data <- cbind(current_exons, intron_ends)
    intron_data$intron_number <- paste("intron", exon_index, sep = "")
    inner_introns <- rbind(inner_introns, intron_data)
  }
  all_introns <- rbind(last_introns, inner_introns)
  all_introns$type <- "intron"

  ## Collecting intron data
  if (verbose) message("Collecting intron data...")
  positive_strand_introns <- subset(all_introns, all_introns$strand == "+")
  negative_strand_introns <- subset(all_introns, all_introns$strand == "-")
  temp_start <- negative_strand_introns$intron_start
  negative_strand_introns$intron_start <- negative_strand_introns$intron_end
  negative_strand_introns$intron_end <- temp_start
  combined_introns <- rbind(positive_strand_introns, negative_strand_introns)
  intron_coordinates <- subset(combined_introns, select = c("seqnames", "intron_start", "intron_end", "width"))
  intron_coordinates$width <- intron_coordinates$intron_end - intron_coordinates$intron_start + 1
  intron_metadata <- combined_introns[, -which(names(combined_introns) %in% c(
    "seqnames", "start", "end", "width", "exon_number", "exon_id",
    "EXON_CLASSIFICATION", "intron_start", "intron_end"
  ))]
  intron_final_data <- cbind(intron_coordinates, intron_metadata)
  grouped_introns <- intron_final_data %>%
    dplyr::group_by(gene_id) %>%
    dplyr::arrange(seqnames, intron_start)
  total_introns <- nrow(grouped_introns)
  if (verbose) message("Total introns: ", total_introns)
  return(as.data.frame(grouped_introns))
}


#' @title Identify Single-Exon Genes/Transcripts in GENCODE Data
#'
#' @description This function identifies single-exon genes or transcripts from a genomic annotation dataset, typically obtained from the GENCODE database. It processes input data, either from a file or an already loaded data frame, and returns either single-exon genes or single-exon transcripts based on the user's selection. The function works seamlessly with GTF/GFF files or data frames that include the standard GENCODE annotation fields. This is especially useful for identifying genes or transcripts that do not undergo splicing.
#'
#' @usage
#' extract_single_exon(input, level = "gene", output_file = NULL)
#'
#'
#' @param input Either a file path to a GTF/GFF3 file or a data frame representing genomic annotations. The input data should contain at least the columns \code{type}, \code{gene_id}, \code{transcript_id}, and \code{exon_number}. If a file path is provided, it is assumed to be a GTF/GFF file.
#'
#' @param level A character string specifying the level of analysis. It should be either \code{"gene"} or \code{"transcript"}. The default is "gene".
#' The selection allows users to specify whether they want to identify single-exon genes or single-exon transcripts in the annotation.
#'
#' @param output_file Optional. A file path to save the results as a tab-separated file. If not provided, the results will not be saved. The file will include the IDs of the single-exon genes or transcripts along with their respective exon IDs. This feature is useful for exporting results for further analysis or reporting purposes.
#'
#' @return A data frame containing the IDs of single-exon genes or transcripts based on the selected level.
#' The returned data frame contains the following columns:
#' - `gene_id` or `transcript_id`: The IDs of the single-exon genes or transcripts.
#' - `exon_count`: The count of exons for each entity (always 1 for single-exon entities).
#' - `exon_id`: The ID of the exon associated with the single-exon gene or transcript.
#'
#' The data frame can then be used for downstream analysis, such as identifying unique single-exon genes or  transcripts, or for exporting to a file.
#'
#' @details
#' \enumerate{
#'   \item Input Validation: Checks for required columns and valid \code{level} specification.
#'   \item Exon Counting: Groups exons by gene/transcript and counts occurrences.
#'   \item Single-Exon Filtering: Retains only entities with exactly one exon.
#'   \item Output: Returns filtered results; optionally writes to file.
#' }
#'
#' @examples
#' # Example input data frame
#' input <- data.frame(
#'   type = c("exon", "exon", "exon", "exon", "exon"),
#'   gene_id = c("gene1", "gene1", "gene2", "gene3", "gene4"),
#'   transcript_id = c("tx1", "tx1", "tx2", "tx3", "tx4"),
#'   exon_number = c(1, 2, 1, 1, 1),
#'   exon_id = c("e1", "e2", "e1", "e1", "e1")
#' )
#'
#' # Identify single-exon genes
#' single_exon_genes <- extract_single_exon(input, level = "gene")
#' print(single_exon_genes)
#'
#' # Identify single-exon transcripts
#' single_exon_transcripts <- extract_single_exon(input, level = "transcript")
#' print(single_exon_transcripts)
#'
#' @seealso \code{\link{extract_introns}}, \code{\link{load_file}}
#'
#' @importFrom utils write.table
#' @export

extract_single_exon <- function(input, level = "gene", output_file = NULL) {
  if (!level %in% c("gene", "transcript")) {
    stop("Invalid level. Please use 'gene' or 'transcript'.")
  }

  # Load input data
  if (is.character(input)) {
    if (!file.exists(input)) {
      stop("The specified file does not exist.")
    }
    message("Reading input file: ", input)
    input_data <- load_file(input)
  } else if (is.data.frame(input)) {
    input_data <- input
  } else {
    stop("Input must be either a file path or a data frame.")
  }
  required_columns <- c("type", "gene_id", "transcript_id", "exon_number", "exon_id")
  if (!all(required_columns %in% colnames(input_data))) {
    stop("Input data must contain the following columns: ", paste(required_columns, collapse = ", "))
  }
  exons <- subset(input_data, input_data$type == "exon")
  id_column <- if (level == "gene") "gene_id" else "transcript_id"
  exon_counts <- as.data.frame(table(exons[[id_column]]))
  colnames(exon_counts) <- c(id_column, "exon_count")

  # Identify single-exon entities
  single_exon <- subset(exon_counts, exon_counts$exon_count == 1)
  single_exon <- merge(single_exon, exons[, c(id_column, "exon_id")], by = id_column, all.x = TRUE)

  # Save results to file if specified
  if (!is.null(output_file)) {
    write.table(single_exon, file = output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    message("Results saved to ", output_file)
  }
  return(single_exon)
}



#' @title Assign intron donor and acceptor splice sites consensus
#'
#' @description This function takes a data frame of intron coordinates and a genome sequence (ideally human or mouse) and returns a data frame with two additional columns for the donor and acceptor splice site consensus sequences.
#' It prepares the donor and acceptor sequences based on the provided intron coordinates and the specified genome (e.g., human hg38), making it useful for downstream analysis of splicing events.
#'
#' @usage
#' assign_splice_sites(input, genome = BSgenome.Hsapiens.UCSC.hg38, verbose = TRUE)
#'
#' @param input A data frame containing intron coordinates with the following columns:
#' - `seqnames`: The chromosome name.
#' - `intron_start`: The start position of the intron.
#' - `intron_end`: The end position of the intron.
#' - `strand`: The strand on which the intron is located (`+` or `-`).
#' - `transcript_id`: The ID of the transcript to which the intron belongs.
#' - `intron_number`: The number of the intron within the transcript.
#' - `gene_name`: The name of the gene.
#' - `gene_id`: The gene ID.
#'
#' @param genome The genome sequence (BSgenome object) for the species. Default is the human genome (hg38).
#' This object is required for extracting the consensus sequences from the genome at the specified intron positions.
#'
#' @param verbose Logical. If TRUE, the function prints progress messages while preparing the splice site data. Default is TRUE.
#'
#'
#' @return A data frame containing the original intron data, with two additional columns:
#' - `donor_ss`: The donor splice site consensus sequence for each intron.
#' - `acceptor_ss`: The acceptor splice site consensus sequence for each intron.
#'
#' @details
#' This function performs the following steps:
#' - First, it prepares the splice site sequences for both donor and acceptor sites by calculating their positions based on the strand orientation and intron coordinates. The donor splice site is typically located at the 5' end of the intron, while the acceptor splice site is at the 3' end.
#' - The function utilizes the `getSeq` function from the `BSgenome` package to extract the nucleotide sequences for both donor and acceptor sites from the specified genome (default is hg38 for humans).
#' - The resulting sequences are added as new columns (`donor_ss` and `acceptor_ss`) to the original input data frame.
#' - The final data frame includes the splice site sequences for each intron, allowing for analysis of splicing efficiency or identification of consensus motifs.
#'
#'
#' @examples
#' suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' gtf_v1 <- load_file(file_v1)
#' introns_df <- extract_introns(gtf_v1)
#' result <- assign_splice_sites(introns_df, genome = BSgenome.Hsapiens.UCSC.hg38)
#'
#' @seealso \code{\link{extract_introns}}, \code{\link{find_cryptic_splice_sites}}
#'
#' @import dplyr
#' @importFrom BSgenome getSeq
#' @export

assign_splice_sites <- function(input, genome = BSgenome.Hsapiens.UCSC.hg38, verbose = TRUE) {
  log_message <- function(msg) {
    if (verbose) cat(paste0("\033[0;32m", msg, "\033[0m\n"))
  }

  required_columns <- c(
    "seqnames", "intron_start", "intron_end", "strand",
    "transcript_id", "intron_number", "gene_name", "gene_id"
  )

  if (!is.data.frame(input) || !all(required_columns %in% colnames(input))) {
    stop(
      "Input must be a data frame containing the following columns: ",
      paste(required_columns, collapse = ", "),
      ". Did you run extract_introns on your input data?"
    )
  }

  # Check if genome is loaded and valid
  if (missing(genome) || is.null(genome)) {
    stop("Genome object is not provided or could not be found.")
  }

  if (!inherits(genome, "BSgenome")) {
    stop("The provided genome object is not a valid BSgenome object.")
  }

  prepare_splice_sites <- function(input, genome, site_type) {
    if (site_type == "donor") {
      input$start_ss <- ifelse(input$strand == "+", input$intron_start, input$intron_end - 1)
      input$end_ss <- ifelse(input$strand == "+", input$intron_start + 1, input$intron_end)
    } else if (site_type == "acceptor") {
      input$start_ss <- ifelse(input$strand == "+", input$intron_end - 1, input$intron_start)
      input$end_ss <- ifelse(input$strand == "+", input$intron_end, input$intron_start + 1)
    }

    # Attempt sequence extraction
    ss_seq <- tryCatch(
      {
        BSgenome::getSeq(
          genome,
          names = input$seqnames,
          start = input$start_ss,
          end = input$end_ss,
          strand = input$strand
        )
      },
      error = function(e) {
        msg <- paste("Failed to retrieve", site_type, "splice site sequences:", e$message)
        stop(msg)
      }
    )

    input[[paste0(site_type, "_ss")]] <- as.character(ss_seq)
    return(input)
  }

  log_message("Preparing donor splice sites data...")
  donor_data <- prepare_splice_sites(input, genome, "donor")

  log_message("Preparing acceptor splice sites data...")
  acceptor_data <- prepare_splice_sites(input, genome, "acceptor")

  log_message("Merging donor and acceptor splice site data...")
  final_data <- input %>%
    dplyr::left_join(donor_data %>% dplyr::select(transcript_id, intron_number, donor_ss),
      by = c("transcript_id", "intron_number")
    ) %>%
    dplyr::left_join(acceptor_data %>% dplyr::select(transcript_id, intron_number, acceptor_ss),
      by = c("transcript_id", "intron_number")
    ) %>%
    dplyr::select(
      seqnames, intron_start, intron_end, strand, transcript_id,
      intron_number, gene_name, gene_id, donor_ss, acceptor_ss
    )

  return(final_data)
}


#' @title Identify Potential Cryptic Splice Sites.
#'
#' @description This function identifies potential cryptic splice sites by comparing sequence motifs in introns to canonical splice site motifs (donor and acceptor). Cryptic splice sites are those that do not match the canonical donor (GT) or acceptor motifs (AG). It compares the identified splice sites with the provided canonical motifs and flags the sites that differ from the canonical patterns, making it useful for studying aberrant splicing events.
#'
#' @usage
#' find_cryptic_splice_sites(input, genome, canonical_donor, canonical_acceptor, verbose)
#'
#'
#' @param input A data frame containing intron coordinates, ideally generated
#'   by `extract_introns()` and `assign_splice_sites()`. Must contain columns: `seqnames`, `intron_start`,
#'   `intron_end`, `strand`, `transcript_id`, `intron_number`, `gene_name`,
#'   `gene_id`, `donor_ss` and `acceptor_ss`.
#' @param genome A BSgenome object representing the genome sequence.
#'   This is used to extract the sequence for each intron to identify splice sites.
#' @param canonical_donor A character vector of canonical donor splice site motifs.
#'   Default is `c("GT")`.
#' @param canonical_acceptor A character vector of canonical acceptor splice site motifs.
#'   Default is `c("AG")`.
#' @param verbose Logical; if `TRUE`, progress messages are printed. Default is `TRUE`.
#'
#' @return The input data frame with two logical columns:
#' \itemize{
#'   \item \code{cryptic_donor}: \code{TRUE} if donor site is non-canonical.
#'   \item \code{cryptic_acceptor}: \code{TRUE} if acceptor site is non-canonical.
#' }
#'
#' @details
#' This function performs the following steps:
#' - It assigns donor and acceptor splice sites to each intron using the `assign_splice_sites` function.
#' - It compares the identified donor and acceptor splice sites against the provided canonical motifs (`GT` for donor and `AG` for acceptor by default). If the splice site sequences do not match the canonical motifs, they are flagged as cryptic.
#' - The function returns a data frame with the same intron information, including additional columns `cryptic_donor` and `cryptic_acceptor` indicating whether the splice sites are cryptic.
#' - The progress of the function is printed if the `verbose` argument is set to `TRUE`, showing also the total number of cryptic donor and acceptor sites and their respective percentages.
#'
#' @examples
#' suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' gtf_v1 <- load_file(file_v1)
#' introns_df <- extract_introns(gtf_v1)
#' introns_ss <- assign_splice_sites(introns_df, genome = BSgenome.Hsapiens.UCSC.hg38)
#' cryptic_sites <- find_cryptic_splice_sites(introns_ss, BSgenome.Hsapiens.UCSC.hg38)
#' head(cryptic_sites)
#'
#' @seealso \code{\link{assign_splice_sites}}, \code{\link{extract_ss_motif}}
#'
#' @importFrom BSgenome getSeq
#' @importFrom dplyr left_join %>% mutate
#' @export

find_cryptic_splice_sites <- function(input, genome,
                                      canonical_donor = c("GT"),
                                      canonical_acceptor = c("AG"),
                                      verbose = TRUE) {
  if (!methods::is(genome, "BSgenome")) {
    stop("The 'genome' argument must be a BSgenome object, but received: ", class(genome))
  }
  if (!is.data.frame(input)) {
    stop("Input must be a data frame. Provided input is of type: ", class(input))
  }
  required_cols <- c(
    "seqnames", "intron_start", "intron_end", "strand",
    "transcript_id", "intron_number", "gene_name", "gene_id",
    "donor_ss", "acceptor_ss"
  )
  missing_cols <- setdiff(required_cols, colnames(input))
  if (length(missing_cols) > 0) {
    stop("Input data frame is missing the following required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (verbose) message("Identifying cryptic splice sites...")
  cryptic_sites <- input %>%
    mutate(
      cryptic_donor = !donor_ss %in% canonical_donor,
      cryptic_acceptor = !acceptor_ss %in% canonical_acceptor
    ) %>%
    filter(cryptic_donor | cryptic_acceptor)
  if (verbose) {
    total_donors <- sum(cryptic_sites$cryptic_donor)
    total_acceptors <- sum(cryptic_sites$cryptic_acceptor)
    total_introns <- nrow(input)
    percent_donors <- (total_donors / total_introns) * 100
    percent_acceptors <- (total_acceptors / total_introns) * 100
    message(sprintf("Detected %d cryptic donors (%.2f%% of total).", total_donors, percent_donors))
    message(sprintf("Detected %d cryptic acceptors (%.2f%% of total).", total_acceptors, percent_acceptors))
  }
  return(cryptic_sites)
}


#' @title Extract Splice Site Motifs for MaxEntScan Analysis (5' or 3')
#'
#' @description This function extracts splice site motifs (5' splice site (5ss) or 3' splice site (3ss)) from a genomic dataset. It retrieves the donor or acceptor splice site motifs for each intron, based on the strand orientation,
#' and compiles them into a FASTA file, which can be used for further analysis (e.g., MaxEntScan).
#'
#' @usage
#' extract_ss_motif(input, genome, type, verbose, save_fasta, output_file)
#'
#' @param input A data frame containing genomic information with the following required columns:
#'   - `seqnames`: Chromosome or scaffold names.
#'   - `strand`: Strand orientation, either '+' or '-'.
#'   - `intron_start`: Start position of the intron.
#'   - `intron_end`: End position of the intron.
#'   - `transcript_id`: Identifier for the transcript.
#'   - `intron_number`: Identifier for the intron within the transcript.
#' @param genome A genome object from the BSgenome package (default is `BSgenome.Hsapiens.UCSC.hg38`).
#' @param type A string indicating which splice site motif to extract. One of `"5ss"` (donor splice site) or `"3ss"` (acceptor splice site).
#' @param verbose Logical; if `TRUE`, progress messages will be printed. Default is `TRUE`.
#' @param save_fasta Logical; if `TRUE`, a FASTA file will be saved containing the extracted motifs. Default is `FALSE`.
#' @param output_file A string specifying the output file path and name for the FASTA file. If `NULL`, a default name will be used (either "5ss_motif_fasta.fa" or "3ss_motif_fasta.fa").
#'
#' @return A data frame with:
#' \itemize{
#'   \item \code{donor_ss_motif} or \code{acceptor_ss_motif}: 9bp (5' ss) or 23bp (3' ss) sequence.
#'   \item Genomic coordinates and transcript metadata.
#' }
#'
#' @details
#' This function performs the following steps:
#' - Based on the `type` argument, the function prepares coordinates for extracting either donor (5ss) or acceptor (3ss) splice site motifs,
#'   adjusting the motif start and end positions depending on the strand orientation.
#' - The motif sequences are then extracted from the specified genome using the `getSeq` function from the BSgenome package.
#' - If `save_fasta` is `TRUE`, a FASTA file is generated containing the extracted motifs, with transcript IDs and intron numbers
#'   used as FASTA headers.
#'
#' @examples
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' gtf_v1 <- load_file(file_v1)
#' introns <- extract_introns(gtf_v1)
#' suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
#' # Extract donor splice site motifs
#' motifs_df <- extract_ss_motif(introns, BSgenome.Hsapiens.UCSC.hg38, "5ss", verbose = FALSE)
#'
#' # Extract acceptor splice site motifs without saving the FASTA file
#' motifs_df <- extract_ss_motif(introns, BSgenome.Hsapiens.UCSC.hg38, "3ss", verbose = FALSE)
#'
#' @seealso \code{\link{assign_splice_sites}}, \code{\link{df_to_fasta}}
#'
#' @importFrom BSgenome getSeq
#' @importFrom dplyr bind_rows
#' @export

extract_ss_motif <- function(input, genome = BSgenome.Hsapiens.UCSC.hg38, type = "5ss", verbose = TRUE, save_fasta = FALSE, output_file = NULL) {
  if (!is.data.frame(input)) {
    stop("'input' must be a dataframe containing genomic data.")
  }
  required_columns <- c("seqnames", "strand", "intron_start", "intron_end", "transcript_id", "intron_number")
  missing_columns <- setdiff(required_columns, colnames(input))
  if (length(missing_columns) > 0) {
    stop(paste("Missing required columns in 'input':", paste(missing_columns, collapse = ", ")))
  }
  if (!type %in% c("5ss", "3ss")) {
    stop("'type' must be either '5ss' (donor splice site) or '3ss' (acceptor splice site).")
  }
  if (type == "5ss") {
    intron_plus <- subset(input, strand == "+")
    intron_plus$motif_start <- intron_plus$intron_start - 3
    intron_plus$motif_end <- intron_plus$intron_start + 5
    intron_minus <- subset(input, strand == "-")
    intron_minus$motif_start <- intron_minus$intron_end - 5
    intron_minus$motif_end <- intron_minus$intron_end + 3
    motif_label <- "donor_ss_motif"
  } else {
    intron_plus <- subset(input, strand == "+")
    intron_plus$motif_start <- intron_plus$intron_end - 19
    intron_plus$motif_end <- intron_plus$intron_end + 3
    intron_minus <- subset(input, strand == "-")
    intron_minus$motif_start <- intron_minus$intron_start - 3
    intron_minus$motif_end <- intron_minus$intron_start + 19
    motif_label <- "acceptor_ss_motif"
  }

  df_plus <- subset(intron_plus, select = c("seqnames", "motif_start", "motif_end", "strand", "transcript_id", "intron_number"))
  df_minus <- subset(intron_minus, select = c("seqnames", "motif_start", "motif_end", "strand", "transcript_id", "intron_number"))
  df <- dplyr::bind_rows(df_plus, df_minus)
  if (verbose) {
    cat(paste0("\033[0;32mPreparing ", motif_label, " sequences ... \033[0m\n"))
  }
  motifs <- as.data.frame(BSgenome::getSeq(genome,
    df$seqnames,
    start = df$motif_start,
    end = df$motif_end,
    strand = df$strand
  ))
  colnames(motifs)[1] <- motif_label
  df2 <- cbind(df, motifs)
  if (save_fasta) {
    if (verbose) {
      cat("Compiling sequences into a FASTA file ...\n")
    }
    if (is.null(output_file)) {
      output_file <- ifelse(type == "5ss", "5ss_motif_fasta.fa", "3ss_motif_fasta.fa")
    }
    df2$fasta_id <- paste0(df2$transcript_id, ";", df2$intron_number)
    df_to_fasta(df2, id_col = "fasta_id", seq_col = motif_label, output_file = output_file, verbose = TRUE)
    df2$fasta_id <- NULL
  }
  return(df2)
}
