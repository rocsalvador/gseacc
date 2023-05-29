#' gseacc
#'
#' A C++ singel sample GSEA for Rna-seq and scRna-seq datasets.
#' As main features:
#' - Mutithreading
#' - Low RAM usage
#'
#' Imports
#' @useDynLib gseacc, .registration = TRUE
#' @importFrom Rcpp
#' @export Gsea
"_PACKAGE"

Rcpp::loadModule(module = "GseaModule", TRUE)
