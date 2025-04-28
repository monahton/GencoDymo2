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
#' \dontrun{
#' # From a GTF file
#' introns <- extract_introns("gencode.v42.gtf")
#'
#' # From a pre-loaded data frame
#' gtf_data <- load_file("gencode.v42.gtf")
#' introns <- extract_introns(gtf_data)
#' }
#'
#' @seealso \code{\link{load_file}}, \code{\link{classify_exons}}
#' @import dplyr
#' @export
#' @keywords GTF introns coordinates splicing

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
