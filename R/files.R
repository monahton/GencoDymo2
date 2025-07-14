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
#' \dontrun{
#' # Retrieve the latest release version for human
#' human_release <- get_latest_release(species = "human", verbose = FALSE)
#' cat("Latest human GENCODE release: release_47")
#'
#' # Get the latest release for mouse
#' mouse_release <- get_latest_release("mouse", verbose = TRUE)
#' }
#'
#' @importFrom RCurl getURL
#' @export

get_latest_release <- function(species, verbose = TRUE) {
  if (!species %in% c("human", "mouse")) {
    stop("Invalid species. Please use 'human' or 'mouse'.")
  }
  base_url <- switch(species,
    human = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
    mouse = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/"
  )
  url <- paste0(base_url, "/")
  ftp_contents <- RCurl::getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
  # Remove HTML tags
  lines <- strsplit(ftp_contents, "\n")[[1]]
  releases <- grep("release_", lines, value = TRUE)
  releases <- sub('.*href="([^"]+)".*', "\\1", releases)
  releases <- trimws(releases)
  releases <- gsub("/", "", releases)
  releases <- releases[grepl("^release_", releases)]
  if (length(releases) == 0) {
    stop("No release directories found at the GENCODE FTP site.")
  }
  latest_release <- releases[length(releases)]
  if (verbose) {
    message(sprintf("Latest %s GENCODE release: %s", species, latest_release))
  }
  return(latest_release)
}

#' @title Download GTF File from the GENCODE Database
#'
#' @description
#'  Downloads a GTF file for a specified species, release version, and annotation type from the GENCODE database.
#' The file is saved to a user-specified directory or the current working directory by default.
#'
#' @usage
#' get_gtf(species, release_version, annotation_type, dest_folder)
#'
#' @param species A character string indicating the species. Supported values are:
#'   - "human"
#'   - "mouse"
#' @param release_version A character string specifying the release version. Options include:
#'   - "latest_release": Automatically fetches the latest release for the specified species.
#'   - "release_X": Specific human release (e.g., "release_47").
#'   - "release_MX": Specific mouse release (e.g., "release_M36").
#' @param annotation_type A character string specifying the annotation type. Valid options are:
#'   - "annotation.gtf.gz"
#'   - "basic.annotation.gtf.gz"
#'   - "chr_patch_hapl_scaff.annotation.gtf.gz"
#'   - "chr_patch_hapl_scaff.basic.annotation.gtf.gz"
#'   - "long_noncoding_RNAs.gtf.gz"
#'   - "primary_assembly.annotation.gtf.gz"
#'   - "primary_assembly.basic.annotation.gtf.gz"
#'   - "tRNAs.gtf.gz"
#'   - "polyAs.gtf.gz"
#' @param dest_folder A character string specifying the destination folder where the file will be downloaded. Defaults to the current working directory.
#'
#' @return A character string specifying the path to the downloaded GTF file.
#'
#' @details
#' The function dynamically determines the correct file URL based on the provided parameters and downloads the GTF file to the desired location.
#' If "latest_release" is specified for `release_version`, the function will first determine the latest available release using `get_latest_release()`.
#'
#'
#' @examples
#' \dontrun{
#' # Download the latest human GTF file with primary assembly annotations into a temp directory
#' temp_dir <- tempdir()
#' gtf_file <- get_gtf(
#'   species = "human",
#'   release_version = "latest_release",
#'   annotation_type = "primary_assembly.basic.annotation.gtf.gz",
#'   dest_folder = temp_dir
#' )
#' print(gtf_file)
#'
#' # Download a specific mouse release with long noncoding RNA annotations into a temp directory
#' temp_dir <- tempdir()
#' gtf_file <- get_gtf(
#'   species = "mouse",
#'   release_version = "release_M36",
#'   annotation_type = "long_noncoding_RNAs.gtf.gz",
#'   dest_folder = temp_dir
#' )
#' print(gtf_file)
#' }
#'
#' @importFrom utils download.file
#' @export

get_gtf <- function(species, release_version = "latest_release", annotation_type, dest_folder = tempdir()) {
  if (!species %in% c("human", "mouse")) {
    stop("Invalid species. Please use 'human' or 'mouse'.")
  }

  if (!release_version %in% c("latest_release") && !grepl("^release_M\\d+$", release_version) && !grepl("^release_\\d+$", release_version)) {
    stop("Invalid release version. Please use 'release_MX' for mouse (e.g., release_M36), or 'release_X' for human, or 'latest_release'.")
  }

  valid_annotation_types <- c(
    "annotation.gtf.gz",
    "basic.annotation.gtf.gz",
    "chr_patch_hapl_scaff.annotation.gtf.gz",
    "chr_patch_hapl_scaff.basic.annotation.gtf.gz",
    "long_noncoding_RNAs.gtf.gz",
    "primary_assembly.annotation.gtf.gz",
    "primary_assembly.basic.annotation.gtf.gz",
    "tRNAs.gtf.gz",
    "polyAs.gtf.gz"
  )

  if (!annotation_type %in% valid_annotation_types) {
    stop("Invalid annotation type. Please use one of the following: ", paste(valid_annotation_types, collapse = ", "))
  }
  base_url <- switch(species,
    human = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
    mouse = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/"
  )

  # Resolve the release version
  if (release_version == "latest_release") {
    release_version <- get_latest_release(species)
  }

  if (species == "mouse") {
    # Mouse uses the vM<number> naming scheme
    version_prefix <- paste0("vM", gsub("release_M", "", release_version))
  } else {
    version_prefix <- paste0("v", gsub("release_", "", release_version))
  }
  file_name <- paste0("gencode.", version_prefix, ".", annotation_type)
  file_url <- paste0(base_url, release_version, "/", file_name)
  local_file <- file.path(dest_folder, file_name)
  message("Downloading GTF file from ", file_url)
  download_result <- tryCatch(
    {
      download.file(file_url, destfile = local_file, mode = "wb")
      TRUE
    },
    error = function(e) {
      message("Error downloading file: ", e$message)
      FALSE
    }
  )
  if (!download_result) {
    stop("Failed to download the GTF file. Please check the URL or your internet connection.")
  }
  message("GTF file downloaded successfully: ", local_file)
  return(local_file)
}


#' @title Download GFF3 File from the GENCODE Database
#'
#' @description
#' Downloads a GFF3 file for a specified species, release version, and annotation type from the GENCODE database.
#' The file is saved to a user-specified directory or the current working directory by default.
#'
#' @usage
#' get_gff3(species, release_version, annotation_type, dest_folder)
#'
#' @param species A character string indicating the species. Supported values are:
#'   - "human"
#'   - "mouse"
#'
#' @param release_version A character string specifying the release version. Options include:
#'   - "latest_release": Fetches the latest release for the species.
#'   - "release_X": Specific release version for human (e.g., "release_42").
#'   - "release_MX": Specific release version for mouse (e.g., "release_M36").
#'
#' @param annotation_type A character string specifying the annotation type. Supported values include:
#'   - "annotation.gff3.gz"
#'   - "basic.annotation.gff3.gz"
#'   - "chr_patch_hapl_scaff.annotation.gff3.gz"
#'   - "chr_patch_hapl_scaff.basic.annotation.gff3.gz"
#'   - "long_noncoding_RNAs.gff3.gz"
#'   - "primary_assembly.annotation.gff3.gz"
#'   - "primary_assembly.basic.annotation.gff3.gz"
#'   - "tRNAs.gff3.gz"
#'   - "polyAs.gff3.gz"
#'
#' @param dest_folder A character string specifying the destination folder. Defaults to the current working directory.
#'
#' @return
#' A character string specifying the full path of the downloaded GFF3 file.
#'
#' @details
#' The function dynamically determines the correct file URL based on the provided parameters and downloads the GFF3 file to the desired location.
#' If "latest_release" is specified for `release_version`, the function will first determine the latest available release using `get_latest_release()`.
#'
#' @examples
#' \dontrun{
#' # Download the latest human GTF file with primary assembly annotations into a temp directory
#' temp_dir <- tempdir()
#' gff3_file <- get_gff3(
#'   species = "human",
#'   release_version = "latest_release",
#'   annotation_type = "primary_assembly.basic.annotation.gff3.gz",
#'   dest_folder = temp_dir
#' )
#' print(gff3_file)
#'
#' # Download a specific mouse release with long noncoding RNA annotations into a temp directory
#' temp_dir <- tempdir()
#' gff3_file <- get_gff3(
#'   species = "mouse",
#'   release_version = "release_M36",
#'   annotation_type = "long_noncoding_RNAs.gff3.gz",
#'   dest_folder = temp_dir
#' )
#' print(gff3_file)
#' }
#'
#' @importFrom utils download.file
#' @export

get_gff3 <- function(species, release_version = "latest_release", annotation_type, dest_folder = tempdir()) {
  if (!species %in% c("human", "mouse")) {
    stop("Invalid species. Please use 'human' or 'mouse'.")
  }

  if (!release_version %in% c("latest_release") && !grepl("^release_M\\d+$", release_version) && !grepl("^release_\\d+$", release_version)) {
    stop("Invalid release version. Please use 'release_MX' for mouse (e.g., release_M36), or 'release_X' for human, or 'latest_release'.")
  }

  valid_annotation_types <- c(
    "annotation.gff3.gz",
    "basic.annotation.gff3.gz",
    "chr_patch_hapl_scaff.annotation.gff3.gz",
    "chr_patch_hapl_scaff.basic.annotation.gff3.gz",
    "long_noncoding_RNAs.gff3.gz",
    "primary_assembly.annotation.gff3.gz",
    "primary_assembly.basic.annotation.gff3.gz",
    "tRNAs.gff3.gz",
    "polyAs.gff3.gz"
  )

  if (!annotation_type %in% valid_annotation_types) {
    stop("Invalid annotation type. Please use one of the following: ", paste(valid_annotation_types, collapse = ", "))
  }
  base_url <- switch(species,
    human = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
    mouse = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/"
  )

  # Resolve the release version
  if (release_version == "latest_release") {
    release_version <- get_latest_release(species)
  }
  if (species == "mouse") {
    # Mouse uses the vM<number> naming scheme
    version_prefix <- paste0("vM", gsub("release_M", "", release_version))
  } else {
    version_prefix <- paste0("v", gsub("release_", "", release_version))
  }
  file_name <- paste0("gencode.", version_prefix, ".", annotation_type)
  file_url <- paste0(base_url, release_version, "/", file_name)
  local_file <- file.path(dest_folder, file_name)
  message("Downloading GFF3 file from ", file_url)
  download_result <- tryCatch(
    {
      download.file(file_url, destfile = local_file, mode = "wb")
      TRUE
    },
    error = function(e) {
      message("Error downloading file: ", e$message)
      FALSE
    }
  )
  if (!download_result) {
    stop("Failed to download the GFF3 file. Please check the URL or your internet connection.")
  }

  message("GFF3 file downloaded successfully: ", local_file)
  return(local_file)
}

#' @title Load a GTF or GFF3 file from GENCODE as a data frame.
#'
#' @description This function imports a GTF or GFF3 file (commonly from the GENCODE website) and converts it into a data frame.
#' The function provides flexibility for users to work with genomic feature files easily in the R environment.
#'
#' @usage load_file(filename)
#' @param filename A character string representing the path to the GTF or GFF3 file (e.g., "gencode.vM36.annotation.gtf.gz").
#' The file could be in GTF or GFF3 format and must be downloaded from a reliable source like GENCODE.
#'
#' @return A data frame containing the parsed content of the GTF or GFF3 file.
#' The data frame includes standard columns such as 'seqnames', 'start', 'end', 'strand', 'feature', and 'gene_id', among others.
#'
#' @details The function uses the `rtracklayer` package to import the GTF or GFF3 file and returns it as a data frame.
#' The user should ensure that the input file is properly formatted and accessible from the specified path.
#' Files larger than a few hundred MBs may take longer to load and process.
#'
#' @examples
#' # Load example GTF files from the package
#' file_v1 <- system.file("extdata", "gencode.v1.example.gtf.gz", package = "GencoDymo2")
#' gtf_v1 <- load_file(file_v1)
#' head(gtf_v1)
#'
#' @importFrom rtracklayer import
#' @importFrom methods is
#' @export

load_file <- function(filename) {
  if (!is.character(filename) || length(filename) != 1) {
    stop("The 'filename' parameter must be a single character string representing the file path.")
  }
  if (!file.exists(filename)) {
    stop("The specified file does not exist. Please check the file path.")
  }
  gtf <- tryCatch(
    {
      rtracklayer::import(filename)
    },
    error = function(e) {
      stop("Error importing the file. Please ensure the file is in GTF or GFF3 format: ", e$message)
    }
  )
  if (!methods::is(gtf, "GRanges")) {
    stop("The imported file is not in a valid GRanges format. Please check the file format.")
  }
  df <- as.data.frame(gtf)
  return(df)
}
