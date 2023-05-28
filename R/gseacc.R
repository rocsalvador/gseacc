#' gseacc
#'
#' A cool package
#'
#' Imports
#' @useDynLib gseacc, .registration = TRUE
#' @importFrom Rcpp
#' @export GseaRcpp
"_PACKAGE"

Rcpp::loadModule(module = "GseaModule", TRUE)
