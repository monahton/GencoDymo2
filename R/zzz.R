.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    paste0("Welcome to GencoDymo2 (v", utils::packageVersion("GencoDymo2"), ")!\n",
           "Tools for dynamic GENCODE annotation analysis.")
  )
}
