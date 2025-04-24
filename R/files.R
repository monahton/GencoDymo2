#' @title Get the Latest Gencode Release Dynamically
#'
#' @description
#' This function retrieves the latest available release of the Gencode database
#' for a given species (human or mouse) by querying the relevant FTP directory.
#' It automates the process of identifying the latest release of annotations for human or mouse.
#'
#' @usage
#' get_latest_release(species, verbose = TRUE)
#'
#' @param species A character string indicating the species. Supported values are:
#'   - "human"
#'   - "mouse"
#' @param verbose Logical. If TRUE (default), the function prints the latest release.
#'
#' @return
#' A character string representing the latest release version for the specified
#' species (e.g., "release_42" for human or "release_36" for mouse).
#'
#' @details
#' The function accesses the GENCODE FTP directory and parses the available folders to determine the latest release.
#' It is particularly useful for bioinformatics workflows that require up-to-date annotation files without manual checks.
#'
#' @examples
#' # Retrieve the latest release version for human
#' human_release <- get_latest_release(species = "human", verbose = FALSE)
#' cat("Latest human GENCODE release: release_47")
#'
#' # Get the latest release for mouse
#' mouse_release <- get_latest_release("mouse", verbose = TRUE)
#'
#' @importFrom RCurl getURL
#' @export
#' @keywords GENCODE, latest release, annotations

get_latest_release <- function(species, verbose = TRUE) {
  if (!species %in% c("human", "mouse")) {
    stop("Invalid species. Please use 'human' or 'mouse'.")
  }
  base_url <- switch(species,
                     human = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
                     mouse = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/")
  url <- paste0(base_url, "/")
  ftp_contents <- RCurl::getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  # Remove HTML tags
  lines <- strsplit(ftp_contents, "\n")[[1]]
  releases <- grep("release_", lines, value = TRUE)
  releases <- sub('.*href="([^"]+)".*', '\\1', releases)
  releases <- trimws(releases)
  releases <- gsub("/", "", releases)
  releases <- releases[grepl("^release_", releases)]
  if (length(releases) == 0) {
    stop("No release directories found at the GENCODE FTP site.")
  }
  release_nums <- as.numeric(gsub("release_", "", releases))
  latest_release <- releases[which.max(release_nums)]
  if (verbose) {
    message(sprintf("Latest %s GENCODE release: %s", species, latest_release))
  }
  return(latest_release)
}




