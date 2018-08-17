chop_slash = function(path) {
  if (stringr::str_sub(path, start=-1) == "/") {
    mod_path = stringr::str_sub(path,end=nchar(path)-1)
  } else {
    mod_path = path
  }
  return(mod_path)
} 

#' Get cytoband data from UCSC
#'
#' @param genome The reference genome (e.g. \code{"hg38"} or \code{"hg19"}).  Default is \code{"hg38"}.
#' @param chr_order The chromosomes of interest and order in which to output the chromosomes.  Default is \code{paste0("chr", c(1:22, "X", "Y"))}.
#'
#' @return A tibble with cytoband data
#' @export
#'
#' @examples
#' cytoband = get_cytobands("hg19", c("chr1", "chr2"))
#' cytoband = get_cytobands() # hg38 and chromosomes 1-22, X and Y.
#' 
get_cytobands <- function(genome = "hg38",
                          chr_order = paste0("chr", c(1:22, "X", "Y"))) {
  genome = "hg38"
  message("opening browser session...")
  session = rtracklayer::browserSession()
  message("creating query...")
  query = rtracklayer::ucscTableQuery(session, "cytoBand", genome)
  message("getting cytoband table...")
  cytoband = rtracklayer::getTable(query)
  colnames(cytoband) = c("chr", "start", "end", "name", "stain")
  cytoband = cytoband %>%
    dplyr::mutate(chr = as.character(chr),
                  start = as.integer(start),
                  end = as.integer(end),
                  name = as.character(name),
                  stain = as.character(stain)) %>%
    dplyr::filter(chr %in% chr_order,
                  name != "") %>%
    dplyr::mutate(chr = factor(chr, levels = chr_order)) %>%
    dplyr::arrange(chr, start) %>%
    tibble::as_tibble()
  message("returning cytoband as tibble...")
  return(cytoband)
}

#' Read in cytoband data (from cytoBandIdeo.txt file downloaded from UCSC)
#'
#' @param cytoband_file cytoBandIdeo.txt file
#' @param chr_order The chromosomes of interest and order in which to output the chromosomes.  Default is \code{paste0("chr", c(1:22, "X", "Y"))}.
#'
#' @return A tibble with cytoband data
#' @export
#'
#' @examples
#' cytoband = read_cytobands("path/to/cytoBandIdeo.txt", c("chr1", "chr2"))
#' cytoband = read_cytobands("path/to/cytoBandIdeo.txt") # chromosomes 1-22, X and Y.
#' 
read_cytobands <- function(cytoband_file,
                           chr_order = paste0("chr", c(1:22, "X", "Y"))) {
  cytoband = readr::read_tsv(cytoband_file,
                             col_names = c("chr", "start", "end", "name", "stain")) %>%
    dplyr::filter(!is.na(name),
                  name != "",
                  chr %in% chr_order) %>%
    dplyr::mutate(chr = factor(chr, levels = chr_order)) %>%
    dplyr::arrange(chr, start)
  return(cytoband)
}

create_chr_mat <- function(all_ordered_data) {
  chrs <- unique(all_ordered_data$loci$chr)
  chr_mat <- sapply(chrs, function(chr) {
    pos = which(all_ordered_data$loci$chr %in% chr)-1
    return(c(pos[1], tail(pos, n=1)))
  })
  chr_mat <- t(chr_mat)
  return(chr_mat)
}

create_chr_arm_mat <- function(loci_info) {
  chr_arms <- unique(loci_info$chr_arm)
  chr_mat <- sapply(chr_arms, function(chr_arm) {
    pos = which(loci_info$chr_arm %in% chr_arm)-1
    return(c(pos[1], tail(pos, n=1)))
  })
  chr_mat <- t(chr_mat)
  return(chr_mat)
}

get_chr_arms = function(cytoband) {
  
  cytoband = cytoband %>%
    dplyr::mutate(arm=substr(name, 1, 1))
  
  chr_arms = dplyr::data_frame(chr=character(),
                               arm=character(),
                               start=integer(),
                               end=integer())
  
  for (i in unique(cytoband$chr)) {
    chr_cytoband = cytoband %>%
      dplyr::filter(chr==i)
    for (j in unique(chr_cytoband$arm)) {
      chr_arm_cytoband = chr_cytoband %>%
        dplyr::filter(arm==j)
      row = dplyr::data_frame(chr=i,
                              arm=j,
                              start=chr_arm_cytoband$start[1],
                              end=chr_arm_cytoband$end[nrow(chr_arm_cytoband)])
      chr_arms = chr_arms %>%
        dplyr::bind_rows(row)
    }
  }
  return(chr_arms)
}

kmeans_init <- function(counts, means, max_states) {
  RCN = counts / (means * sum(counts))
  out = Ckmeans.1d.dp::Ckmeans.1d.dp(RCN, k = max_states)
  return(list(state_rcns=out$centers,
              hidden_states=out$cluster-1))
}

init_hidden_states <- function(counts, means, max_states, kmeans=TRUE) {
  if(kmeans) {
    return(kmeans_init(counts, means, max_states))
  } else {
    return(list(state_rcns=rep(0, max_states),
                hidden_states=rep(0, length(counts))))
  }
}

order_data <- function(loci, counts, chr_order=paste0("chr", c(1:22, "X", "Y"))) {
  loci$chr <- factor(loci$chr, levels=chr_order)
  loci_order <- with(loci,order(chr,leftPos))
  loci_ordered <- loci[loci_order,]
  loci_ordered$chr <- droplevels(loci_ordered$chr)
  counts_ordered <- counts[loci_order,]
  loci_counts <- dplyr::bind_cols(loci_ordered, counts_ordered)
  return(list(counts=counts_ordered,
              loci=loci_ordered,
              loci_counts=loci_counts))
}
