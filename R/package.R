#' Genealogical Simulations for Two-Segment Viruses
#' 
#' @name SimReassort-package
#' @docType package
#' @author Qianying Lin
#' 
#' @import ggplot2
#' @importFrom foreach registerDoSEQ
#' 
#' @useDynLib phylopomp, .registration = TRUE, .fixes="P_"
#' 
NULL

foreach::registerDoSEQ()
