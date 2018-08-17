getH = function(CN_vec, gridsize) {
  out <- tryCatch(
    {
      KernSmooth::dpik(CN_vec, gridsize=gridsize)
    },
    error=function(cond) {
      return(KernSmooth::dpik(CN_vec, gridsize=gridsize, scalest = "stdev"))
    },
    finally={}
  )
  return(out)
}

getMAP = function(CN_vec) {
  if (length(unique(CN_vec)) == 1) {
    return(CN_vec[1])
  } else {
    gridsize = length(CN_vec)
    #h = dpik(CN_vec, gridsize = gridsize)
    h = getH(CN_vec, gridsize)
    dens = KernSmooth::bkde(CN_vec, bandwidth=h, gridsize=gridsize)
    return(MAP = dens$x[which.max(dens$y)])
  }
}

calculate_populated_states <- function(hidden_states) {
  sapply(hidden_states, function(x) {length(unique(x))})
}

Mode <- function(x) {
  # Note this only returns the first Mode.
  # For example, for the vector: c(1,2) it will return 1.
  y <- unique(x)
  return(y[which.max(tabulate(match(x, y)))])
}

# Functions to read in MCMC chains

read_MCMC_run <- function(path, sample_name) {
  path = chop_slash(path)
  sample_path <- paste0(path, "/", sample_name)
  sample1 <- list()
  sample1$params <- readr::read_tsv(paste0(sample_path,"_params.txt"),
                                    col_names=TRUE)
  sample1$stateRCN <- readr::read_tsv(paste0(sample_path,"_cn.txt"),
                                      col_names=FALSE,
                                      col_types = readr::cols(.default=readr::col_double()))
  sample1$hs <- readr::read_tsv(paste0(sample_path,"_hs.txt"),
                                col_names=FALSE,
                                col_types = readr::cols(.default=readr::col_character()))
  sample1$precision <- readr::read_tsv(paste0(sample_path,"_p.txt"),
                                       col_names="prec", 
                                       col_types="d")
  sample1$alpha <- readr::read_tsv(paste0(sample_path,"_alpha.txt"),
                                   col_names="alpha",
                                   col_types="d")
  sample1$gamma <- readr::read_tsv(paste0(sample_path,"_gamma.txt"),
                                   col_names="gamma",
                                   col_types="d")
  sample1$kappa <- readr::read_tsv(paste0(sample_path,"_kappa.txt"), 
                                   col_names="kappa",
                                   col_types="d")
  sample1$mcll <- readr::read_tsv(paste0(sample_path,"_mcll.txt"),
                                  col_names="mcll",
                                  col_types="d")
  sample1$ibbsll <- readr::read_tsv(paste0(sample_path,"_ibbsll.txt"),
                                    col_names="ibbsll", 
                                    col_types="d")
  sample1$beta <- readr::read_tsv(paste0(sample_path,"_beta.txt"), 
                                  col_names=FALSE, 
                                  col_types = readr::cols(.default=readr::col_double()))
  return(sample1)
}

read_MCMC_hs <- function(path, sample_name) {
  path = chop_slash(path)
  sample_path <- paste0(path, "/", sample_name)
  hs <- readr::read_tsv(paste0(sample_path,"_hs.txt"),
                        col_names=FALSE,
                        col_types = readr::cols(.default=read::col_character()))
  return(hs)
}

read_MCMC_params <- function(path, sample_name) {
  path = chop_slash(path)
  sample_path <- paste0(path, "/", sample_name)
  sample_info <- readr::read_tsv(paste0(sample_path,"_params.txt"),
                                 col_names=TRUE)
  return(sample_info)
}

read_MCMC_HP_and_precision <- function(path, sample_name) {
  sample_path <- paste0(path, sample_name)
  sample1 <- list()
  sample1$params <- readr::read_tsv(paste0(sample_path,"_params.txt"),
                                    col_names=TRUE)
  sample1$precision <- readr::read_tsv(paste0(sample_path,"_p.txt"),
                                       col_names="prec",
                                       col_types="d")
  sample1$alpha <- readr::read_tsv(paste0(sample_path,"_alpha.txt"),
                                   col_names="alpha", 
                                   col_types="d")
  sample1$gamma <- readr::read_tsv(paste0(sample_path,"_gamma.txt"), 
                                   col_names="gamma", 
                                   col_types="d")
  sample1$kappa <- readr::read_tsv(paste0(sample_path,"_kappa.txt"), 
                                   col_names="kappa", 
                                   col_types="d")
  sample1$mcll <- readr::read_tsv(paste0(sample_path,"_mcll.txt"),
                                  col_names="mcll", 
                                  col_types="d")
  sample1$ibbsll <- readr::read_tsv(paste0(sample_path,"_ibbsll.txt"), 
                                    col_names="ibbsll", 
                                    col_types="d")
  return(sample1)
}
