#' FAST-SeqS aligned read counts for OAC, GAC, BO and control samples
#' 
#' A dataset containing the FAST-SeqS aligned read counts for oesophageal adenocarcinoma, 
#' Barrett's oesophagus, gastric adenocarcinoma and control samples.  These are
#' the counts obtained by the pipeline described in Abujudeh et al. 2018. Low-cost and clinically applicable
#' copy number profiling using repeat DNA. bioRxiv doi: 
#'
#' @format A data frame with 215 columns (of which 212 are samples) and 58725 rows representing genomic loci
#' 
#' \describe{
#'   \item{chr}{Chromosome}
#'   \item{leftPos}{The genomic coordinate within `chr`}
#'   \item{strand}{The strand, '+' or '-'}
#' }
#' @references Abujudeh, S. et al. (2018) Low-cost and clinically applicable
#' copy number profiling using repeat DNA. \emph{bioRxiv doi: }
"loci_counts"