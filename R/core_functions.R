
#' Choose and filter loci with zeros
#' 
#' Given a list of samples, counts data frame 
#' (with columns for \code{chr}, \code{leftPos} and \code{strand}),
#' and chromosomes of interest
#' return a data frame with loci filtered 
#' (loci with less than or equal to 
#' \code{zeros_permitted} in the counts of the 
#' \code{samples_used_to_filter} are kept).
#'
#' @param analysis_path The path to the directory in which 
#' you want the results to be saved.
#' If the path does not exist, it will be created for you.
#' @param counts A data frame with \code{chr},
#' \code{leftPos} and \code{strand} 
#' in the first three columns and sample counts in all other columns
#' @param chr_order A vector of chromosomes of interest.
#' e.g. \code{paste0("chr", 1:22)} or 
#' \code{paste0("chr", c(1:22, "X", "Y"))}
#' @param samples_used_to_filter A character vector with sample names
#'  used to filter the counts
#' @param zeros_permitted Loci with less than or equal to 
#' \code{zeros_permitted} in the counts of the 
#' \code{samples_used_to_filter} are kept.  Default: 0.
#' @param output_folder The folder name to store 
#' the filtered counts (Default: `filtered_counts`).
#' Not recommended to change this.
#' @param counts_Rdata The name of the 
#' resulting Rdata object (Default: \code{loci_counts.Rdata}). 
#' Not recommended to change this.
#'
#' @return A counts data frame (after filtering loci with 
#' zero counts in any of the samples provided)
#' @export
#'
#' @examples
#' filtered_counts = choose_loci("/path/to/analysis",
#' loci_counts,
#' c("sample_1", "sample_2", "sample_3"))
#' 
#' # see data(loci_counts) for example counts data frame
choose_loci <- function(analysis_path,
                        counts,
                        chr_order,
                        samples_used_to_filter,
                        zeros_permitted = 0,
                        output_folder="filtered_counts",
                        counts_Rdata="loci_counts.Rdata") {

  if (missing(analysis_path)) {
    analysis_path = getwd()
  }
  
  if (stringr::str_sub(analysis_path, start=-1) == "/") {
    analysis_path = stringr::str_sub(analysis_path,
                            end=nchar(analysis_path)-1)
  }
  
  if (!dir.exists(analysis_path)) {
    message("analysis_path directory doesn't exist, creating directory: ", 
            file.path(analysis_path, output_folder))
    dir.create(analysis_path)
  }
  
  if (!dir.exists(file.path(analysis_path, output_folder))) {
    message("creating directory to store output: ", 
            file.path(analysis_path, output_folder))
    dir.create(file.path(analysis_path, output_folder))
  }
  

  if (!all.equal(names(counts)[1:3], c("chr", "leftPos", "strand"))) {
    stop("Must provide counts data frame with the first 3 columns: 'chr', 'leftPos', 'strand'")
  }
  

  message("Keeping the following chromosomes:\n", 
          paste(shQuote(chr_order), collapse="\n"))
  
  counts = counts %>%
    dplyr::filter(chr %in% chr_order) %>%
    dplyr::mutate(chr = factor(chr, levels = chr_order)) %>%
    dplyr::arrange(chr, leftPos)
    
  message("Using the following samples to filter loci:\n", 
          paste(shQuote(samples_used_to_filter), collapse="\n"))
  sample_counts = counts %>%
    dplyr::select(dplyr::one_of(samples_used_to_filter))
  
  message("Filtering loci with less than or equal to ",
          zeros_permitted, " zeros in any of the samples...")
  loci_robust = apply(sample_counts, 1, function(y) {
    sum(y==0) <= zeros_permitted
  })
  
  counts = counts %>%
    dplyr::filter(loci_robust)
  
  message("Saving filtered counts in Rdata object: ",
          file.path(analysis_path, output_folder, counts_Rdata))
  save(counts, file=file.path(analysis_path, output_folder, counts_Rdata))
  message("Returning counts...")
  return(counts)
}

#' Fit the expected proportion of reads at each locus using a set of control samples
#' 
#' This is usually run after you have filtered the loci with zero counts 
#' (using \code{\link{choose_loci}} function).  It is important that, for each locus, 
#' at least one control sample has a non-zero count.
#'
#' @param analysis_path The path to the directory in which 
#' you want the results to be saved.  If the path does
#'  not exist, it will be created for you.
#' @param counts A data frame of counts with \code{chr},
#' \code{leftPos} and \code{strand}
#' in the first three columns and sample counts in all other columns
#' @param sample_types A data frame with two columns
#' The first column should have the title \code{Sample} and contain 
#' the name of all the samples you wish to use in this analysis 
#' (this is typically all the muliplexed samples on a sequencing lane).
#' The second column should have the title \code{Control} and the 
#' value of each row should be 1 if the sample is a control 
#' sample or 0 if it is not.
#' @param output_folder The folder name to store the fit data
#' (Default: `fit_controls`). 
#' This will be created for you if it doesn't exist already.
#' Not recommended to change this.
#' @param output_run_name The name to prefix the output files.
#' Default: `fit_controls`. Not recommended to change this.
#' @param iterations The number of iterations to run the MCMC.
#' Default: \code{20000}.
#' @param burn_in The number of initial iterations to burn. 
#' Default: \code{5000}.
#' @param precision_prior The Gamma parameters for the prior 
#' on the sample inverse dispersion parameter, parameterised by 
#' shape and scale (Default: \code{c(1.5, 1000000)}).
#' @param mean_priors The Beta prior parameters for the proportion
#'  of reads at a locus in a control sample.
#' This is represented by a matrix with 2 rows and L columns 
#' (where L is the number of rows in your counts file).
#' The first row represents the alpha shape parameter and 
#' the second row represents the beta shape parameter.
#' Default: \code{matrix(c(1,1), nrow=2, ncol=nrow(counts))}
#' i.e. a flat `Beta(1,1)` prior for each locus. 
#' @param precision_sigma The standard deviation of the proposal 
#' jump distribution for the sample inverse dispersion parameter. 
#' Default: \code{50000}.
#' @param mean_sigma The standard deviation of the proposal jump
#'  distribution for the sample inverse dispersion parameter. 
#'  Default: \code{0.0001}.
#'
#' @return The output is a list with 1) a vector with the mean 
#' of the proportion values (\code{mean_means}), 2) MAP estimate of 
#' the proportion values (\code{map_means}), 3) A list with count and 
#' loci data for the samples of interest (\code{run_data})
#' @export
#'
#' @examples
#' sample_types = data.frame(Samples=c("sample_1", "sample_2", "sample_3"), Control=c(1, 1, 0)) # in practice you'll want more samples.
#' fit_data = fit_controls(analysis_path = "/path/to/analysis", filtered_counts, sample_types)
fit_controls <- function(analysis_path,
                         counts,
                         sample_types,
                         output_folder,
                         output_run_name,
                         iterations,
                         burn_in,
                         precision_prior, # shape and scale of gamma prior
                         mean_priors,
                         precision_sigma,
                         mean_sigma) {

  if (missing(analysis_path)) {
    analysis_path = getwd()
  }
  
  if (stringr::str_sub(analysis_path, start=-1) == "/") {
    analysis_path = stringr::str_sub(analysis_path,
                            end=nchar(analysis_path)-1)
  }
  
  num_loci = nrow(counts)
  
  if (missing(output_folder)) {
    output_folder = "fit_controls"
  }
  
  if (missing(output_run_name)) {
    output_run_name = "fit_controls"
  }
  
  if (missing(iterations)) {
    iterations = 20000
  }
  
  if (missing(burn_in)) {
    burn_in = 5000
  }
  
  if (missing(precision_prior)) {
    precision_prior = c(1.5, 1000000)
  }
  
  if (missing(mean_priors)) {
    mean_priors = matrix(1, ncol=num_loci, nrow=2)
  }
  
  if (missing(precision_sigma)) {
    precision_sigma = 50000
  }
  
  if (missing(mean_sigma)) {
    mean_sigma = 0.0001
  }
  
  if (!dir.exists(analysis_path)) {
    dir.create(analysis_path)
  }
  
  if (!dir.exists(file.path(analysis_path, output_folder))) {
    dir.create(file.path(analysis_path, output_folder))
  }

  control_samples = sample_types %>%
    dplyr::filter(Control == 1) %>%
    .$Sample
  message("control samples are:\n", paste0(control_samples, "\n"))
  
  control_counts = counts %>%
    dplyr::select(dplyr::one_of(control_samples))
  
  control_count_matrix = control_counts %>%
    as.matrix()
  
  all_samples = sample_types %>%
    .$Sample
  message("all samples are:\n", paste0(all_samples, "\n"))

  # store old working directory.
  # This is a temporary measure until we update 
  # inferBBParamsWithPriors to include an output directory
  # for now, we set the output directory as the working directory.
  message("running MCMC to infer loci proportions...")
  old_dir = getwd()
  setwd(file.path(analysis_path,output_folder))
  system.time(inferBBParamsWithPriors(control_count_matrix,
                                      iterations,
                                      mean_sigma, 
                                      precision_sigma,
                                      mean_priors,
                                      precision_prior,
                                      output_run_name))
  
  message("reading in MCMC data....")
  inferred_means <- matrix(scan(paste0(output_run_name,
                                       "_means.txt")),
                           ncol=num_loci, byrow=TRUE)
  
  message("selecting iterations after burn in (>=",
          burn_in,
          " iterations)...")
  inferred_means <- inferred_means[burn_in:iterations,]
  
  message("summarizing loci proportions...")
  mean_means <- colMeans(inferred_means)
  mean_means = mean_means / sum(mean_means)
  
  map_means = apply(inferred_means, 2, getMAP)
  map_means = map_means / sum(map_means)
  
  run_loci = counts %>%
    dplyr::select(chr, leftPos, strand)
  
  counts = counts %>%
    dplyr::select(-chr, -leftPos, -strand)

  run_data = order_data(run_loci, counts)
  save(mean_means, map_means, run_data,
       file=paste0(output_run_name, ".Rdata"))
  file.remove(paste0(output_run_name, "_means.txt"))
  setwd(old_dir)
  
  # TODO: add here typical diagnostic plots of the inference of the means and precision
  # TODO: output an "updated prior" for the precision of new samples
  
  return(list(mean_means=mean_means,
              map_means=map_means,
              run_data=run_data))
}


#' Run Markov Chain Monte Carlo (MCMC) to infer RCN profiles
#'
#' This function runs the MCMC algorithm and outputs files
#' with samples of the posterior distribution for:
#' loci hidden states, state RCN, 
#' sample precision/inverse dispersion parameter
#' and other latent variables in the conliga model. This function should
#' be run after \code{\link{choose_loci}} and \code{\link{fit_controls}}.
#' Use \code{\link{process_MCMC}} to process the results of
#' \code{\link{run_MCMC}}.  
#'
#' @param analysis_path The path to the directory in 
#' which you want the results to be saved.
#' If the path does not exist, it will be created for you.
#' @param fit_data The fit_data list created by the \code{\link{fit_controls}} function
#' @param cytoband Cytoband information provided by the 
#' \code{\link{get_cytobands}} or \code{\link{read_cytobands}} functions
#' @param priors A data frame with the profiles defined.
#' If not provided, defaults will be used.
#' @param samples A vector of sample names you 
#' wish to infer RCN profiles for
#' @param suffix Optional: if you want to add a suffix
#'  to the sample name
#' @param num_cores The number of cores to use. Default: 1.
#' If more than one core is used, the \code{mcmapply} function
#' from the \code{parallel} package
#' is used.
#' @param output_folder The output folder for the MCMC chains.
#' Default is `MCMC`.
#' @param iterations The number of iterations to the run the MCMC for
#' (Default: \code{50000})
#' @param thin The thinning of the chain
#' (Default: \code{10})
#' @param gamma The hyperparameter gamma (must be > 0).
#' Default: \code{1}.
#' If you prefer to infer gamma from the data 
#' (not recommended!) set it to \code{-1}.
#' @param rho The hyperparameter rho (value between 0-1). 
#' This controls how "sticky" the hidden Markov chain is.
#' Sensible values are around \code{0.99}. Default: \code{-1}, 
#' meaning it is inferred from the data.
#' @param alpha_plus_kappa The hyperparameter 
#' alpha_plus_kappa (must be > 0).
#' This controls the variability of the transition matrix
#' The greater the value, the less variance between 
#' rows of the transition matrix
#' Default: \code{-1}, meaning its value is inferred from the data
#' @param max_states The value to use to truncate 
#' the Dirichlet Process.
#' The default is \code{30}.  Do not set this too low.
#' It must be greater than the number of expected 
#' copy number states,
#' otherwise the Dirichlet approximation to the Dirichlet 
#' Process will not hold. 
#' @param init_states Boolean: \code{TRUE}/\code{FALSE}.
#' If \code{TRUE}, 1d k-means clustering will be performed 
#' prior to running the chain (using the \code{Ckmeans.1d.dp}
#' package). This ensures the chain starts in a reasonable place.
#' If \code{FALSE}, the chain will be started in a random starting point.
#' Note by choosing \code{FALSE}, there is greater risk 
#' that the MCMC will fail in the first few iterations 
#' of the MCMC due to starting in a very unlikely position.
#' The default is \code{TRUE}.
#'
#' @return Returns \code{NULL} and creates MCMC output files 
#' for each sample in \code{analysis_path/output_folder/}.
#' @export
#'
#' @examples
#' run_MCMC(analysis_path=analysis_path,
#' fit_data=fit_data,
#' cytoband=cytoband,
#' priors=priors,
#' samples=samples_to_run,
#' num_cores=num_cores,
#' iterations=50000,
#' thin=10,
#' gamma = 1, # define gamma set to 1,
#' max_states=30,
#' init_states=TRUE) # start chain off by k-means clustering
#' 
#' # see \link{https://github.com/samabs/conliga} for more information
run_MCMC <- function(analysis_path, 
                     fit_data,
                     cytoband,
                     priors,
                     samples,
                     suffix,
                     num_cores,
                     output_folder,
                     iterations,
                     thin,
                     gamma,
                     rho,
                     alpha_plus_kappa,
                     max_states,
                     init_states) {

  if (missing(analysis_path)) {
    analysis_path = getwd()
  }
  
  if (stringr::str_sub(analysis_path, start=-1) == "/") {
    analysis_path = stringr::str_sub(analysis_path,
                            end=nchar(analysis_path)-1)
  }
  
  if (missing(output_folder)) {
    output_folder = "MCMC"
  }
  
  if (missing(iterations)) {
    iterations = 50000
  }
  
  if (missing(thin)) {
    thin = 10
  }
  
  if (missing(suffix)) {
    suffix = rep("", length(samples))
  }
  
  if (missing(num_cores)) {
    num_cores = 1
  }
  
  if (!dir.exists(analysis_path)) {
    dir.create(analysis_path)
  }
  
  if (!dir.exists(file.path(analysis_path, 
                            output_folder))) {
    dir.create(file.path(analysis_path, 
                         output_folder))
  }
  
  if (missing(gamma)) {
    gamma = 1 # invalid value and treated by MCMC as gamma not given
  } else {
    if(is.numeric(gamma) & !(gamma > 0)) {
      stop("gamma must be greater than 0")
    } else {
      message("gamma provided as fixed with value: ",
              gamma)
    }
  }
  
  if (missing(rho)) {
    rho = -1 # invalid value and treated by MCMC as rho not given
  } else {
    if(is.numeric(rho) & (rho < 0 | rho > 1)) {
      stop("rho must be between 0 and 1")
    } else {
      message("rho provided as fixed with value: ",
              rho)
    }
  }
  
  if (missing(alpha_plus_kappa)) {
    alpha_plus_kappa = -1 # invalid value and treated by MCMC as rho not given
  } else {
    if(is.numeric(alpha_plus_kappa) & !(alpha_plus_kappa > 0)) {
      stop("alpha_plus_kappa must be greater than 0")
    } else {
      message("alpha_plus_kappa provided as fixed with value: ",
              alpha_plus_kappa)
    }
  }
  
  if (missing(max_states)) {
    max_states = 30
  }
  
  if (missing(init_states)) {
    init_states = rep(TRUE, length(samples))
  }
  
  if (length(init_states)==1) {
    init_states = rep(init_states, length(samples))
  }
  
  chr_arms = get_chr_arms(cytoband)
  
  chr_arm_ranges <- GenomicRanges::GRanges(seqnames=chr_arms$chr,
                                           ranges=IRanges::IRanges(start=chr_arms$start,
                                                                   end=chr_arms$end),
                                           arm=chr_arms$arm)
  
  loci_ranges <- GenomicRanges::GRanges(seqnames=fit_data$run_data$loci$chr,
                                        ranges=IRanges::IRanges(start=fit_data$run_data$loci$leftPos,
                                                                end=fit_data$run_data$loci$leftPos+1))
  
  overlaps <- GenomicRanges::findOverlaps(loci_ranges, chr_arm_ranges)
  
  loci_arm <- rep(NA, times=length(nrow(fit_data$run_data$loci)))
  loci_arm[GenomicRanges::queryHits(overlaps)] <- chr_arms$arm[GenomicRanges::subjectHits(overlaps)]
  
  loci_arms <- fit_data$run_data$loci
  loci_arms <- loci_arms %>%
    dplyr::mutate(arm=loci_arm)
  loci_arms <- loci_arms %>%
    dplyr::mutate(chr_arm=paste0(loci_arms$chr, loci_arms$arm))
  
  chr_arm_mat <- create_chr_arm_mat(loci_arms)
  
  message("processing the following samples: \n", 
          paste0(samples, suffix, "\n"))

  
  if (missing(priors)) {
    priors = dplyr::data_frame(prior=c("hpp_rho_c",
                                       "hpp_rho_d",
                                       "hpp_g_a",
                                       "hpp_g_b",
                                       "hpp_apk_a",
                                       "hpp_apk_b",
                                       "precision_prior_shape",
                                       "precision_prior_scale",
                                       "rcn_prior_shape",
                                       "rcn_prior_scale"),
                               value=c(100000,
                                       100,
                                       7,
                                       1,
                                       2000,
                                       10,
                                       1.5,
                                       1000000,
                                       3,
                                       1))
  } 
  
  if (num_cores > 1) {
    num_samples = length(samples)
    
    # make variables to run MCMC in parallel
    hpp_rho_c = rep(priors %>%
                      dplyr::filter(prior == "hpp_rho_c") %>%
                      .[["value"]], num_samples)
    hpp_rho_d = rep(priors %>%
                      dplyr::filter(prior == "hpp_rho_d") %>%
                      .[["value"]], num_samples)
    hpp_g_a = rep(priors %>%
                    dplyr::filter(prior == "hpp_g_a") %>%
                    .[["value"]], num_samples)
    hpp_g_b = rep(priors %>%
                    dplyr::filter(prior == "hpp_g_b") %>%
                    .[["value"]], num_samples)
    hpp_apk_a = rep(priors %>%
                      dplyr::filter(prior == "hpp_apk_a") %>%
                      .[["value"]], num_samples)
    hpp_apk_b = rep(priors %>%
                      dplyr::filter(prior == "hpp_apk_b") %>%
                      .[["value"]], num_samples)
    precision_prior_shape = rep(priors %>%
                                  dplyr::filter(prior == "precision_prior_shape") %>%
                                  .[["value"]], num_samples)
    precision_prior_scale = rep(priors %>%
                                  dplyr::filter(prior == "precision_prior_scale") %>%
                                  .[["value"]], num_samples)
    rcn_prior_shape = rep(priors %>%
                            dplyr::filter(prior == "rcn_prior_shape") %>%
                            .[["value"]], num_samples)
    rcn_prior_scale = rep(priors %>%
                            dplyr::filter(prior == "rcn_prior_scale") %>%
                            .[["value"]], num_samples)
    
    gamma_list = rep(list(gamma), num_samples)
    rho_list = rep(list(rho), num_samples)
    apk_list = rep(list(alpha_plus_kappa), num_samples)
    
    chr_arm_mat_list = rep(list(chr_arm_mat), num_samples)
    means_list = rep(list(fit_data$map_means), num_samples)
    
    counts_list = lapply(fit_data$run_data$counts[samples], unlist)
    
    old_dir = getwd()
    setwd(file.path(analysis_path,output_folder))
    
    init_data = samples %>%
      purrr::map2(.y = 1:length(samples),
                  ~init_hidden_states(fit_data$run_data$counts[[.x]],
                                      fit_data$map_means,
                                      max_states,
                                      init_states[.y])) %>%
      purrr::transpose()
    
    if(all(init_states)) {
      message("starting MCMC in states by k-means clustering...")
    } else {
      message("starting MCMC in randomly sampled states for some samples...")
      message(paste(samples, suffix, " ", init_states, 
                    " (TRUE=k-means, FALSE=random)\n"))
      warning("WARNING: starting MCMC in random start may lead to failure of MCMC in first few iterations!")
    }
    
    message("running MCMC in parallel with ", 
            num_cores, 
            " cores...")
    run_sticky_in_parallel(num_cores,
                           counts_list,
                           chr_arm_mat_list,
                           means_list,
                           hpp_g_a,
                           hpp_g_b,
                           hpp_apk_a,
                           hpp_apk_b,
                           hpp_rho_c,
                           hpp_rho_d,
                           paste0(samples, suffix),
                           precision_prior_shape,
                           precision_prior_scale,
                           rcn_prior_shape,
                           rcn_prior_scale,
                           gamma_list,
                           rho_list,
                           apk_list,
                           init_hidden_states_list = init_data$hidden_states,
                           init_state_rcns_list = init_data$state_rcns,
                           iterations=iterations,
                           thin=thin,
                           max_states=max_states)
    setwd(old_dir)
  } else {

    old_dir = getwd()
    setwd(file.path(analysis_path,output_folder))
    
    if(all(init_states)) {
      message("starting MCMC in states by k-means clustering...")
    } else {
      message("starting MCMC in randomly sampled states for some samples...")
      message(paste(samples, suffix, " ",
                    init_states, " (TRUE=k-means, FALSE=random)\n"))
      warning("WARNING: starting MCMC in random start may lead to failure of MCMC in first few iterations!")
    }
    
    purrr::pmap(list(samples, suffix, 1:length(samples)),
                function(s, p, n) {
      
      init_data = init_hidden_states(fit_data$run_data$counts[[s]],
                                     fit_data$map_means,
                                     max_states,
                                     init_states[n])
      
      run_sticky(counts=fit_data$run_data$counts[[s]],
                 chr_mat=chr_arm_mat,
                 means=fit_data$map_means,
                 iterations=iterations,
                 gamma_a=priors$value[priors$prior=="rcn_prior_shape"],
                 gamma_scale=priors$value[priors$prior=="rcn_prior_scale"],
                 sample_ref=paste0(s, p),
                 thin=thin,
                 hpp_g_a=priors$value[priors$prior=="hpp_g_a"],
                 hpp_g_b=priors$value[priors$prior=="hpp_g_b"],
                 hpp_ak_a=priors$value[priors$prior=="hpp_apk_a"],
                 hpp_ak_b=priors$value[priors$prior=="hpp_apk_b"],
                 hpp_r_c=priors$value[priors$prior=="hpp_rho_c"],
                 hpp_r_d=priors$value[priors$prior=="hpp_rho_d"],
                 precision_p_shape=priors$value[priors$prior=="precision_prior_shape"],
                 precision_p_scale=priors$value[priors$prior=="precision_prior_scale"],
                 gamma=gamma,
                 rho=rho,
                 alpha_plus_kappa=alpha_plus_kappa,
                 max_states=max_states,
                 init_hidden_states=init_data$hidden_states,
                 init_state_rcns=init_data$state_rcns)
    })

    setwd(old_dir)
  }
}

#' Process and summarise the MCMC output from \code{\link{run_MCMC}}
#' 
#' This function process the output from \code{\link{run_MCMC}}.
#' It produces standard plots along with the relative copy number profile 
#' and credible intervals for \code{sample_name}.  It outputs this as a \code{tsv}
#' file in \code{analysis_path/results_output_folder/sample_name/sample_name.tsv}.
#' 
#'
#' @param analysis_path The path to the directory in 
#' which you want the results to be saved.
#' If the path does not exist, it will be created for you.
#' @param sample_name The sample name for which you want to obtain results
#' @param fit_data The \code{fit_data} list produced by \code{\link{fit_controls}}
#' @param cytoband Cytoband data frame produced by 
#' \code{\link{get_cytobands}} or \code{\link{read_cytobands}} functions
#' @param burn_in The first iterations to burn from the MCMC chain (Default: 5000).
#' For example, if the MCMC was run for 50,000 \code{iterations}
#'  with a \code{thin} value of 10, a \code{burn_in} of 5000 would 
#'  discard the first 5000 iterations (500 iterations aftering thinning)
#'  leaving you with 45000 iterations (4500 after thinning)
#' @param MCMC_output_folder The directory relative to \code{analysis_path} 
#' where the MCMC data was output by \code{\link{run_MCMC}}.  
#' The default is `MCMC`` and we recommend not to change this.
#' @param results_output_folder The folder, relative to \code{analysis_path} in
#' which the summarised output and plots will be saved.  Default is "results"
#' and we recommend not to change this.
#'
#' @return Returns \code{NULL} and creates some basic plots 
#' and a \code{sample_name.tsv} file
#' in \code{analysis_path/results_output_folder/sample_name/}.
#' @export
#'
#' @examples
#' process_MCMC("/path/to/analysis/",
#' "sample_1", fit_data, cytoband, 5000)
#' 
#' # see \link{https://github.com/samabs/conliga} for more information
#' 
process_MCMC <- function(analysis_path,
                         sample_name,
                         fit_data,
                         cytoband,
                         burn_in,
                         MCMC_output_folder,
                         results_output_folder) {

  if (missing(analysis_path)) {
    analysis_path = getwd()
  }
  
  if (stringr::str_sub(analysis_path, start=-1) == "/") {
    analysis_path = stringr::str_sub(analysis_path,
                            end=nchar(analysis_path)-1)
  }
  
  if (missing(MCMC_output_folder)) {
    MCMC_output_folder = "MCMC"
  }
  
  if (missing(results_output_folder)) {
    results_output_folder = "results"
  }
  
  if (missing(burn_in)) {
    burn_in = 5000
  }
    
  if (!dir.exists(analysis_path)) {
    dir.create(analysis_path)
  }
  
  if (!dir.exists(file.path(analysis_path,
                            results_output_folder))) {
    dir.create(file.path(analysis_path,
                         results_output_folder))
  }
  
  if (!dir.exists(file.path(analysis_path, 
                            results_output_folder, 
                            sample_name))) {
    dir.create(file.path(analysis_path, 
                         results_output_folder, 
                         sample_name))
  }
  
  if (!dir.exists(file.path(analysis_path, 
                            results_output_folder,
                            sample_name, 
                            "standard_plots"))) {
    dir.create(file.path(analysis_path, 
                         results_output_folder, 
                         sample_name, 
                         "standard_plots"))
  }
  
  message("processing MCMC output for sample: ", 
          sample_name)
  chr_arms = get_chr_arms(cytoband)
  message("created chromosome arms...")
  sample_params <- read_MCMC_params(file.path(analysis_path,MCMC_output_folder), 
                                    sample_name)
  currentwd = getwd()
  setwd(file.path(analysis_path,MCMC_output_folder))
  p_file = paste0(sample_name,"_probs.bin")
  n_states_file = paste0(sample_name,"_nstates.bin")

  message("getting state switching permutations using MAP number of states method...")
  # get permutations for map number of states
  perm_map = get_state_permutations(p_file, n_states_file, sample_params[["iterations"]],
                                    burn_in, sample_params[["thin"]],
                                    sample_params[["max_states"]], sample_params[["num_loci"]],
                                    TRUE, 100)

  message("getting state switching permutations using MAXIMUM number of states method...")
  # get permutations for maximum number of states
  perm_max = get_state_permutations(p_file, n_states_file, sample_params[["iterations"]],
                                    burn_in, sample_params[["thin"]],
                                    sample_params[["max_states"]], sample_params[["num_loci"]],
                                    FALSE, 100)

  # read in MCMC data.
  message("reading in the rest of the MCMC data...")
  MCMC_data = read_MCMC_run(file.path(analysis_path,MCMC_output_folder), sample_name)
  
  # observed iterations of MCMC
  message("iterations of interest...")
  iter_w_data <- base::seq(from=0, to=sample_params[["iterations"]], by=sample_params[["thin"]])
  
  
  message("organising state assignments...")
  # state assignments
  sa = data.frame(MCMC_data$hs,
                  iter = iter_w_data)
  colnames(sa) = c(paste0("l", 1:sample_params[["num_loci"]]), "iter")
  sa = sa %>% tidyr::gather(loci, index, -iter)
  sa = sa %>%
    dplyr::filter(iter > burn_in) %>%
    dplyr::mutate(index = as.character(index))

  message("organising state relative copy number data...")
  # state relative copy number
  srcn = data.frame(MCMC_data$stateRCN,
                    iter = iter_w_data)
  colnames(srcn) = c(0:(sample_params[["max_states"]]-1), "iter")
  srcn = srcn %>%
    tidyr::gather(index, CN, -iter)
  srcn = srcn %>%
    dplyr::filter(iter > burn_in) %>%
    dplyr::mutate(index = as.character(index))
  
  message("joining state assignments and state relative copy number data...")
  srcn = sa %>%
    dplyr::left_join(srcn)
  
  message("Summarising relative copy number by loci...")
  srcn_summarise = srcn %>%
    dplyr::group_by(loci) %>%
    dplyr::summarise(map_rcn = getMAP(CN),
                     min_rcn = min(CN),
                     max_rcn = max(CN),
                     # mean_rcn = mean(CN),
                     q2p5 = quantile(CN, probs=0.025),
                     q5 = quantile(CN, probs=0.05),
                     q10 = quantile(CN, probs=0.1),
                     q20 = quantile(CN, probs=0.2),
                     q30 = quantile(CN, probs=0.3),
                     q40 = quantile(CN, probs=0.4),
                     q50 = quantile(CN, probs=0.5),
                     q60 = quantile(CN, probs=0.6),
                     q70 = quantile(CN, probs=0.7),
                     q80 = quantile(CN, probs=0.8),
                     q90 = quantile(CN, probs=0.9),
                     q95 = quantile(CN, probs=0.95),
                     q97p5 = quantile(CN, probs=0.975)) %>%
    dplyr::ungroup()
  
  chr_order = unique(fit_data$run_data$loci$chr)
  
  message("Joining with loci location information...")
  chr_pos = fit_data$run_data$loci %>%
    dplyr::mutate(chr = factor(chr, levels=chr_order))
  chr_pos = chr_pos %>%
    dplyr::mutate(loci = paste0("l",1:nrow(chr_pos)))
  
  srcn_summarise = chr_pos %>%
    dplyr::right_join(srcn_summarise)
  srcn_summarise = srcn_summarise %>%
    dplyr::mutate(chr = factor(chr, levels=chr_order)) %>%
    dplyr::arrange(chr, leftPos)
  
  message("Getting probability of loci state changes...")
  hs = MCMC_data$hs %>%
    dplyr::mutate(iter = iter_w_data)
  hs = hs %>%
    dplyr::filter(iter > burn_in) %>%
    dplyr::select(-iter) %>%
    as.matrix()
  class(hs) = "numeric" # to convert from character to numeric
  diff = lociStateDiff(hs)
  dim(diff) = NULL # to transform into an R vector
  diff = diff / nrow(hs)
  
  message("Converting coords to Mbp in chromosome arms...")
  chr_arms = chr_arms %>%
    dplyr::mutate(chr = factor(chr, levels=chr_order),
                  start_Mbp=start / 1000000,
                  end_Mbp=end / 1000000) 

  message("Making changepoint data frame...")
  change_point = data.frame(change=diff, 
                            loci=paste0("l",1:length(diff)))
  change_point = change_point %>%
    dplyr::left_join(chr_pos) %>%
    dplyr::mutate(chr = factor(chr, levels=chr_order))
  
  change_point = change_point %>%
    dplyr::group_by(chr) %>%
    dplyr::mutate(change = ifelse(dplyr::row_number() == 1, 0, change)) %>%
    dplyr::ungroup()

  message("Making change point plots...")
  pp = change_point %>% 
    dplyr::mutate(leftPos = leftPos / 1000000) %>%
    ggplot2::ggplot(ggplot2::aes(x=leftPos, y = change)) +
    ggplot2::geom_line() +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=start_Mbp[1]), 
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=end_Mbp), 
                        colour="blue", 
                        linetype="dashed") +
    theme_sticky() +
    ggplot2::scale_x_continuous(labels=scales::comma) +
    ggplot2::xlab("chromsome position (Mbp)") +
    ggplot2::ylab("state change probability") +
    ggplot2::facet_wrap(~chr, scales="free_x")
  
  ggplot2::ggsave(file.path(analysis_path, 
                            results_output_folder, 
                            sample_name, 
                            "standard_plots",
                            paste0(sample_name, 
                                   "_change_point_l.png")),
                  pp, height=15, width=20, dpi = 150)
  pp = change_point %>% 
    dplyr::mutate(leftPos = leftPos / 1000000) %>%
    ggplot2::ggplot(ggplot2::aes(x=leftPos, y = change)) +
    ggplot2::geom_point(stroke=0) +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=start_Mbp[1]), 
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=end_Mbp), 
                        colour="blue", 
                        linetype="dashed") +
    theme_sticky() +
    ggplot2::scale_x_continuous(labels=scales::comma) +
    ggplot2::xlab("chromsome position (Mbp)") +
    ggplot2::ylab("state change probability") +
    ggplot2::facet_wrap(~chr, scales="free_x")
  
  ggplot2::ggsave(file.path(analysis_path, 
                            results_output_folder, 
                            sample_name, 
                            "standard_plots",
                            paste0(sample_name, 
                                   "_change_point_p.png")),
                  pp, height=15, width=20, dpi = 150)
  
  message("making RCN line plots...")
  pp = srcn_summarise %>%
    dplyr::mutate(leftPos = leftPos / 1000000) %>%
    ggplot2::ggplot(ggplot2::aes(x=leftPos, y=map_rcn)) +
    ggplot2::geom_ribbon(ggplot2::aes(x=leftPos, 
                                      ymin=q2p5, 
                                      ymax=q97p5), 
                         fill="#D3D3D3") +
    ggplot2::geom_ribbon(ggplot2::aes(x=leftPos, 
                                      ymin=q5, 
                                      ymax=q95), 
                         fill="#BEBEBE") +
    ggplot2::geom_ribbon(ggplot2::aes(x=leftPos, 
                                      ymin=q10, 
                                      ymax=q90), 
                         fill="#A9A9A9") +
    ggplot2::geom_line() +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=start_Mbp[1]), 
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=end_Mbp), 
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::scale_x_continuous(labels=scales::comma) +
    ggplot2::xlab("chromsome position (Mbp)") +
    ggplot2::ylab("relative copy number") +
    ggplot2::facet_wrap(~chr, scales="free") +
    theme_sticky()
  
  ggplot2::ggsave(file.path(analysis_path, 
                            results_output_folder, 
                            sample_name, 
                            "standard_plots",
                            paste0(sample_name, 
                                   "_rcn_line_scaled.png")),
                  pp, height=15, width=20, dpi = 150)
  
  pp = srcn_summarise %>%
    dplyr::mutate(leftPos = leftPos / 1000000) %>%
    ggplot2::ggplot(ggplot2::aes(x=leftPos, 
                                 y=map_rcn)) +
    ggplot2::geom_ribbon(ggplot2::aes(x=leftPos, 
                                      ymin=q2p5, 
                                      ymax=q97p5), 
                         fill="#D3D3D3") +
    ggplot2::geom_ribbon(ggplot2::aes(x=leftPos, 
                                      ymin=q5, 
                                      ymax=q95), 
                         fill="#BEBEBE") +
    ggplot2::geom_ribbon(ggplot2::aes(x=leftPos, 
                                      ymin=q10, 
                                      ymax=q90), 
                         fill="#A9A9A9") +
    ggplot2::geom_line() +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=start_Mbp[1]),
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=end_Mbp), 
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::scale_x_continuous(labels=scales::comma) +
    ggplot2::xlab("chromsome position (Mbp)") +
    ggplot2::ylab("relative copy number") +
    ggplot2::facet_wrap(~chr, scales="free_x") +
    theme_sticky()
  
  ggplot2::ggsave(file.path(analysis_path, 
                            results_output_folder, 
                            sample_name, 
                            "standard_plots",
                            paste0(sample_name, 
                                   "_rcn_line_unscaled.png")),
                  pp, height=15, width=20, dpi = 150)

  message("making RCN point plots...")
  pp = srcn_summarise %>%
    dplyr::mutate(leftPos = leftPos / 1000000) %>%
    ggplot2::ggplot(ggplot2::aes(x=leftPos, 
                                 y=map_rcn)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=q2p5, 
                                        ymax=q97p5), 
                           alpha=0.6, 
                           colour="grey", 
                           size=0.2, 
                           width=3) +
    ggplot2::geom_point(colour="black", 
                        stroke=0) +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=start_Mbp[1]), 
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=end_Mbp), 
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::scale_x_continuous(labels=scales::comma) +
    ggplot2::xlab("chromsome position (Mbp)") +
    ggplot2::ylab("relative copy number") +
    ggplot2::facet_wrap(~chr, scales="free") +
    theme_sticky()
  ggplot2::ggsave(file.path(analysis_path, 
                            results_output_folder, 
                            sample_name, 
                            "standard_plots",
                            paste0(sample_name, 
                                   "_rcn_point_scaled.png")),
                  pp, height=15, width=20, dpi = 150)
  
  pp = srcn_summarise %>%
    dplyr::mutate(leftPos = leftPos / 1000000) %>%
    ggplot2::ggplot(ggplot2::aes(x=leftPos, 
                                 y=map_rcn)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=q2p5, 
                                        ymax=q97p5), 
                           alpha=0.6, 
                           colour="grey", 
                           size=0.2, 
                           width=3) +
    ggplot2::geom_point(colour="black", 
                        stroke=0) +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=start_Mbp[1]),
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=end_Mbp), 
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::scale_x_continuous(labels=scales::comma) +
    ggplot2::xlab("chromsome position (Mbp)") +
    ggplot2::ylab("relative copy number") +
    ggplot2::facet_wrap(~chr, scales="free_x") +
    theme_sticky()
  
  ggplot2::ggsave(file.path(analysis_path, 
                            results_output_folder, 
                            sample_name, 
                            "standard_plots",
                            paste0(sample_name, 
                                   "_rcn_point_unscaled.png")),
                  pp, height=15, width=20, dpi = 150)

  # read in number of states in each iteration
  message("Reading in number of states in each MCMC iteration...")
  con <- file(paste0(file.path(analysis_path,MCMC_output_folder), 
                     "/", sample_name, "_nstates.bin"), 
              'rb')
  n_states = readBin(con, "integer", 
                     n=sample_params[["iterations"]] / sample_params[["thin"]], 
                     size=4)
  close(con)

  message("keeping only the states after burn in...")  
  burned_n_states = n_states[((burn_in/sample_params[["thin"]])+1):(sample_params[["iterations"]] / sample_params[["thin"]])]
  
  nn = srcn %>%
    dplyr::group_by(iter) %>%
    dplyr::summarise(num_states=dplyr::n_distinct(index))

  num_loci = srcn %>%
    dplyr::group_by(iter, index) %>%
    dplyr::summarise(num_loci = dplyr::n_distinct(loci))

  message("swapping labels for MAP solution...")
  labels_map = swap_labels(MCMC_data, 
                           perm_map, 
                           burned_n_states, 
                           nn, 
                           n_states, 
                           iter_w_data, 
                           burn_in, 
                           TRUE)

  message("swapping labels for MAX solution...")
  labels_max = swap_labels(MCMC_data, 
                           perm_max, 
                           burned_n_states,
                           nn, 
                           n_states, 
                           iter_w_data, 
                           burn_in, 
                           FALSE)
  
  # remove states in iterations where there were fewer than max states
  labels_max = labels_max %>%
    dplyr::filter(!is.na(index))

  message("removing MCMC data from memory...")
  rm(MCMC_data)
  gc()
  
  colnames(labels_map)[2] = "new_state_map"
  colnames(labels_max)[2] = "new_state_max"

  message("making state / iteration plots...")
  pp1 = labels_map %>%
    dplyr::left_join(num_loci) %>%
    ggplot2::ggplot(ggplot2::aes(x=iter, 
                                 y=CN, 
                                 size=num_loci)) +
    ggplot2::geom_point(ggplot2::aes(colour=index), 
                        stroke=0, 
                        alpha=0.5) +
    ggplot2::ylab("RCN") +
    ggplot2::scale_size() +
    theme_sticky()
  
  pp2 = labels_map %>%
    dplyr::left_join(num_loci) %>%
    ggplot2::ggplot(ggplot2::aes(x=iter, 
                                 y=CN, 
                                 size=num_loci)) + 
    ggplot2::geom_point(ggplot2::aes(colour=new_state_map), 
                        stroke=0, 
                        alpha=0.5) +
    ggplot2::ylab("RCN") +
    ggplot2::scale_size() +
    theme_sticky()
  
  pp3 = labels_max %>%
    dplyr::left_join(num_loci) %>%
    ggplot2::ggplot(ggplot2::aes(x=iter,
                                 y=CN, 
                                 size=num_loci)) +
    ggplot2::geom_point(ggplot2::aes(colour=index), 
                        stroke=0, 
                        alpha=0.5) +
    ggplot2::ylab("RCN") +
    ggplot2::scale_size() +
    theme_sticky()
  
  pp4 = labels_max %>%
    dplyr::left_join(num_loci) %>%
    ggplot2::ggplot(ggplot2::aes(x=iter, 
                                 y=CN, 
                                 size=num_loci)) + 
    ggplot2::geom_point(ggplot2::aes(colour=new_state_max), 
                        stroke=0, 
                        alpha=0.5) +
    ggplot2::ylab("RCN") +
    ggplot2::scale_size() +
    theme_sticky()

  message("Made label switching diagnosis plots, now putting them together...")
  pp = gridExtra::arrangeGrob(pp1, pp2, pp3, pp4)
  message("Saving plots...")
  ggplot2::ggsave(file.path(analysis_path, 
                            results_output_folder,
                            sample_name, 
                            "standard_plots",
                            paste0(sample_name, 
                                   "_states_rcn_switching.png")),
                  pp, height=15, width=20, dpi = 150)
  
  message("Making map_map data frame...")

  max_map = labels_max %>%
    dplyr::filter(!is.na(CN))

  max_map = max_map %>%
    dplyr::group_by(new_state_max)

  max_map = max_map %>%
    dplyr::summarise(max_state_CN = getMAP(CN),
                     max_state_CN_q97p5 = quantile(CN, probs=0.975),
                     max_state_CN_q2p5 = quantile(CN, probs=0.025))

 
  message("Making map_map data frame")
  map_map = labels_map %>%
    dplyr::filter(!is.na(CN)) %>%
    dplyr::group_by(new_state_map) %>%
    dplyr::summarise(map_state_CN = getMAP(CN),
                     map_state_CN_q97p5 = quantile(CN, probs=0.975),
                     map_state_CN_q2p5 = quantile(CN, probs=0.025))
  
  message("Making plots with MAP state estimates")
  pp2 = labels_map %>%
    dplyr::left_join(map_map) %>%
    dplyr::group_by(iter, new_state_map) %>%
    dplyr::group_by(iter, new_state_map) %>%
    dplyr::summarise(CN = dplyr::first(CN),
                     map_state_CN = dplyr::first(map_state_CN),
                     map_state_CN_q97p5 = dplyr::first(map_state_CN_q97p5),
                     map_state_CN_q2p5 = dplyr::first(map_state_CN_q2p5)) %>%
    ggplot2::ggplot(ggplot2::aes(x=iter, y=CN, colour=new_state_map)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=map_state_CN, 
                                     colour=new_state_map)) + 
    ggplot2::geom_point(alpha=0.4, stroke=0) +
    theme_sticky()
  
  pp1 = labels_map %>%
    dplyr::group_by(iter, index) %>%
    dplyr::summarise(CN = dplyr::first(CN)) %>%
    ggplot2::ggplot(ggplot2::aes(x=iter, 
                                 y=CN, 
                                 colour=index)) +
    ggplot2::geom_point(alpha=0.4, stroke=0) +
    theme_sticky()
  
  pp3 = labels_max %>%
    dplyr::group_by(iter, index) %>%
    dplyr::summarise(CN = dplyr::first(CN)) %>%
    ggplot2::ggplot(ggplot2::aes(x=iter, 
                                 y=CN, 
                                 colour=index)) +
    ggplot2::geom_point(alpha=0.4, stroke=0) + 
    theme_sticky()
  
  pp4 = labels_max %>%
    dplyr::left_join(max_map) %>%
    dplyr::group_by(iter, new_state_max) %>%
    dplyr::summarise(CN = dplyr::first(CN),
                     max_state_CN = dplyr::first(max_state_CN),
                     max_state_CN_q97p5 = dplyr::first(max_state_CN_q97p5),
                     max_state_CN_q2p5 = dplyr::first(max_state_CN_q2p5)) %>%
    ggplot2::ggplot(ggplot2::aes(x=iter, 
                                 y=CN, 
                                 colour=new_state_max)) +
    ggplot2::geom_point(alpha=0.4, stroke=0) +
    ggplot2::geom_hline(ggplot2::aes(yintercept=max_state_CN, 
                                     colour=new_state_max)) + 
    theme_sticky()
  
  message("Putting plots together...")
  pp = gridExtra::arrangeGrob(pp1, pp2, pp3, pp4)
  
  message("Saving plot...")
  ggplot2::ggsave(file.path(analysis_path, 
                            results_output_folder, 
                            sample_name, 
                            "standard_plots",
                   paste0(sample_name, 
                          "_states_rcn_switching_estimate.png")),
         pp, height=15, width=20, dpi = 150)
  
  
  message("summarising data by new state labels...")
  new_labels = labels_map %>%
    dplyr::select(-CN) %>%
    dplyr::right_join(labels_max)
  new_labels = new_labels %>%
    dplyr::select(-CN) %>%
    dplyr::left_join(srcn)
  
  new_labels_map = new_labels %>%
    dplyr::filter(!is.na(new_state_map)) %>%
    dplyr::group_by(loci) %>%
    dplyr::summarise(new_state_map = Mode(new_state_map)) %>%
    dplyr::left_join(chr_pos) %>%
    dplyr::left_join(map_map) %>%
    dplyr::arrange(chr, leftPos)
  
  new_labels_max = new_labels %>%
    dplyr::filter(!is.na(new_state_max)) %>%
    dplyr::group_by(loci) %>%
    dplyr::summarise(new_state_max = Mode(new_state_max)) %>%
    dplyr::left_join(chr_pos) %>%
    dplyr::left_join(max_map) %>%
    dplyr::arrange(chr, leftPos)
  
  rm(new_labels)
  
  pp = new_labels_map %>%
    dplyr::mutate(leftPos = leftPos / 1000000) %>%
    ggplot2::ggplot(ggplot2::aes(x=leftPos, 
                                 y=map_state_CN)) +
    ggplot2::geom_point(ggplot2::aes(colour=new_state_map), 
                        stroke=0) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=map_state_CN_q97p5, 
                                        ymax=map_state_CN_q2p5, 
                                        colour=new_state_map), 
                           alpha=0.2) +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=start_Mbp[1]), 
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=end_Mbp),
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::scale_x_continuous(labels=scales::comma) +
    ggplot2::xlab("chromsome position (Mbp)") +
    ggplot2::ylab("relative copy number") +
    ggplot2::facet_wrap(~chr, scales="free_x") +
    theme_sticky()
  
  message("saving relative copy number plots WITH uncertainty information and with using new labels (MAP label switch method)...")
  ggplot2::ggsave(file.path(analysis_path, 
                            results_output_folder,
                            sample_name, 
                            "standard_plots",
                            paste0(sample_name, 
                                   "_MAP_solution_error.png")),
                  pp, height=15, width=20, dpi = 150)
  
  pp = new_labels_map %>%
    dplyr::mutate(leftPos = leftPos / 1000000) %>%
    ggplot2::ggplot(ggplot2::aes(x=leftPos, y=map_state_CN)) +
    ggplot2::geom_point(ggplot2::aes(colour=new_state_map), 
                        stroke=0) +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=start_Mbp[1]), 
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=end_Mbp), 
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::scale_x_continuous(labels=scales::comma) +
    ggplot2::xlab("chromsome position (Mbp)") +
    ggplot2::ylab("relative copy number") +
    ggplot2::facet_wrap(~chr, scales="free_x") +
    theme_sticky()
  
  message("saving relative copy number plots WITHOUT uncertainty information and with using new labels (MAP label switch method)...")
  ggplot2::ggsave(file.path(analysis_path, 
                            results_output_folder, 
                            sample_name, 
                            "standard_plots",
                            paste0(sample_name, 
                                   "_MAP_solution.png")),
                  pp, height=15, width=20, dpi = 150)
  
  pp = new_labels_max %>%
    dplyr::mutate(leftPos = leftPos / 1000000) %>%
    dplyr::filter(!is.na(chr)) %>%
    ggplot2::ggplot(ggplot2::aes(x=leftPos, 
                                 y=max_state_CN)) +
    ggplot2::geom_point(ggplot2::aes(colour=new_state_max), 
                        stroke=0) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=max_state_CN_q97p5, 
                                        ymax=max_state_CN_q2p5, 
                                        colour=new_state_max), 
                           alpha=0.2) +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=start_Mbp[1]), 
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=end_Mbp), 
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::scale_x_continuous(labels=scales::comma) +
    ggplot2::xlab("chromsome position (Mbp)") +
    ggplot2::ylab("relative copy number") +
    ggplot2::facet_wrap(~chr, scales="free_x") +
    theme_sticky()
  
  message("saving relative copy number plots WITH uncertainty information and with using new labels (MAX label switch method)...")
  ggplot2::ggsave(file.path(analysis_path, 
                            results_output_folder,
                            sample_name, 
                            "standard_plots",
                   paste0(sample_name,
                          "_MAX_solution_error.png")),
         pp, height=15, width=20, dpi = 150)
  
  pp = new_labels_max %>%
    dplyr::mutate(leftPos = leftPos / 1000000) %>%
    dplyr::filter(!is.na(chr)) %>%
    ggplot2::ggplot(ggplot2::aes(x=leftPos, 
                                 y=max_state_CN)) +
    ggplot2::geom_point(ggplot2::aes(colour=new_state_max),
                        stroke=0) +
    ggplot2::geom_vline(data=chr_arms, 
                        ggplot2::aes(xintercept=start_Mbp[1]),
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::geom_vline(data=chr_arms,
                        ggplot2::aes(xintercept=end_Mbp),
                        colour="blue", 
                        linetype="dashed") +
    ggplot2::scale_x_continuous(labels=scales::comma) +
    ggplot2::xlab("chromsome position (Mbp)") +
    ggplot2::ylab("relative copy number") +
    ggplot2::facet_wrap(~chr, scales="free_x") +
    theme_sticky()
  
  message("saving relative copy number plots WITHOUT uncertainty information and with using new labels (MAX label switch method)...")
  ggplot2::ggsave(file.path(analysis_path,
                            results_output_folder, 
                            sample_name, 
                            "standard_plots",
                            paste0(sample_name, 
                                   "_MAX_solution.png")),
                  pp, height=15, width=20, dpi = 150)
  
  
  all_data = srcn_summarise %>%
    dplyr::left_join(new_labels_map) %>%
    dplyr::left_join(new_labels_max) %>%
    dplyr::left_join(change_point)
  
  message("saving TSV file with summarised results...")
  readr::write_tsv(all_data, 
                   file.path(analysis_path, 
                             results_output_folder, 
                             sample_name,
                             paste0(sample_name, 
                                    "_results.tsv")))
  return(NULL)
}

