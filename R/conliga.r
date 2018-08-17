#' conliga
#' 
#' This package calls relative copy number profiles using FAST-SeqS data (or similar)
#' 
#' @details
#' See the GitHub page (\link{https://github.com/samabs/conliga}) for common functions and examples.
#' 
#' The main functions and workflow:
#' 
#' \enumerate{
#' \item{\code{\link{choose_loci}}}
#' \item{\code{\link{fit_controls}}}
#' \item{\code{\link{run_MCMC}}}
#' \item{\code{\link{process_MCMC}}}
#' }
#' 
#' More documentation to come.
#' 
#' @docType package
#' @author Sam Abujudeh <sam.abujudeh@cruk.cam.ac.uk>
#' @importFrom Rcpp evalCpp
#' @useDynLib conliga
#' @name conliga
"_PACKAGE"