
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
--------

Installation
------------

``` r
# Obtain the development version of conliga from Github
# install.packages("devtools")
devtools::install_github("samabs/conliga")
```

Usage
-----

Load the `conliga` package.

``` r
library(conliga)
```

The typical workflow:

1.  `choose_loci()`: filter loci.
2.  `fit_controls()`: fit a vector of proportions, using a set of control samples
3.  `run_MCMC()`: run the MCMC to obtain posterior probabilities of the latent variables in the conliga model
4.  `process_MCMC()`: process the MCMC to summarise the posterior probabilities to obtain relative copy number profiles and basic plots

Example
-------

Let's go through an example using the data from the paper.

### Set-up

First let's load the loci count data.

``` r
# load the loci count data for the samples
data(loci_counts)

loci_counts[1:6, 1:6]
#>    chr leftPos strand SLX-9394.FastSeqA SLX-9394.FastSeqB
#> 1 chr1  800730      +                 0                 0
#> 2 chr1  882779      -                 3                13
#> 3 chr1 2366307      +                 0                 0
#> 4 chr1 3617459      +                 0                 0
#> 5 chr1 3954290      -                 0                 0
#> 6 chr1 4045043      +                 0                 0
#>   SLX-9394.FastSeqC
#> 1                 0
#> 2                23
#> 3                 0
#> 4                 0
#> 5                 0
#> 6                 0
```

Now we have a `loci_counts` data frame with 58725 rows corresponding to genomic loci. The data frame contains 215 columns. The first 3 columns contain `chr`, `leftPos` and `strand` and correspond to the chromosome, position and strand of a locus. The other 212 columns represent count observations for each sample.

### Filtering loci using `choose_loci()`

Now we filter the rows of the data frame to filter out noise and retain 'robust' loci for further analysis. We choose a set of samples to use as a basis for filtering and a threshold on the number of zeros permitted in the set of samples. It can be useful to list the samples we want to use as a basis for filtering in a csv or tsv file.

``` r
samples_file = system.file("extdata",
                           "samples_for_robust_loci.csv",
                           package = "conliga",
                           mustWork = TRUE)

samples_for_filtering = readr::read_csv(samples_file) %>%
  .$Sample

samples_for_filtering
#>  [1] "SLX-9394.FastSeqE" "SLX-9394.FastSeqF" "SLX-9394.FastSeqG"
#>  [4] "SLX-9394.FastSeqH" "SLX-9394.FastSeqI" "SLX-9394.FastSeqJ"
#>  [7] "SLX-9394.FastSeqK" "SLX-9395.FastSeqA" "SLX-9395.FastSeqB"
#> [10] "SLX-9395.FastSeqC" "SLX-9395.FastSeqD" "SLX-9395.FastSeqE"
#> [13] "SLX-9395.FastSeqF" "SLX-9395.FastSeqG" "SLX-9395.FastSeqH"
```

At the very least, we must ensure that at least one of the control samples used for fitting the proportions (see next section) has a non-zero count. In this example, we aggressively filter the loci by removing a locus with a zero count in any of out control samples.

Now we run the `choose_loci` function, passing it: \* an analysis path which acts as out base directory for all output in our analysis. \* the counts data frame \* the chromosomes we are interested in (in this case chr1-22) \* the samples we use as a basis for filtering \* the number of permitted zeros

If a locus has more than the number of permitted zeros in the samples list and/or it is in a chromosome we are not interested in, it is filtered.

``` r
# the base directory for the analysis output files
analysis_path = "~/my_analysis_path"
 
# the chromosomes we want to analyse
chromosomes = paste0("chr", 1:22)

filtered_counts = choose_loci(analysis_path=analysis_path,
                              counts=loci_counts,
                              chr_order=chromosomes, 
                              samples_used_to_filter=samples_for_filtering,
                              zeros_permitted=0)

filtered_counts[1:6, 1:6]
#> # A tibble: 6 x 6
#>   chr   leftPos strand `SLX-9394.FastSe… `SLX-9394.FastS… `SLX-9394.FastS…
#>   <fct>   <int> <chr>              <int>            <int>            <int>
#> 1 chr1   882779 -                      3               13               23
#> 2 chr1  4045044 +                     11               23               31
#> 3 chr1  4080807 +                      8               16               41
#> 4 chr1  4309374 +                      1                1                1
#> 5 chr1  4350131 +                      3                3               20
#> 6 chr1  4404254 -                     12               14               35
```

Now we are left with 11854 loci in the `filtered_counts` data frame.

Fitting the proportion of read counts using `fit_controls()`
------------------------------------------------------------

Now we use a set of controls (which are assumed to be diploid at each locus) from a batch, to fit the proportion of read counts at each locus.

We need to provide a data frame which describes the samples we will be analysing and defines the samples we will use as controls:

``` r
sample_types_file = system.file("extdata",
                                "samples_info.csv",
                                package = "conliga",
                                mustWork = TRUE)

sample_types = readr::read_csv(sample_types_file)
sample_types
#> # A tibble: 24 x 2
#>    Sample            Control
#>    <chr>               <int>
#>  1 SLX-9394.FastSeqC       0
#>  2 SLX-9394.FastSeqD       0
#>  3 SLX-9394.FastSeqE       1
#>  4 SLX-9394.FastSeqG       1
#>  5 SLX-9394.FastSeqH       1
#>  6 SLX-9394.FastSeqI       1
#>  7 SLX-9394.FastSeqJ       1
#>  8 SLX-9394.FastSeqK       1
#>  9 SLX-9395.FastSeqA       1
#> 10 SLX-9395.FastSeqC       1
#> # ... with 14 more rows
```

The column `Control` is `1` if the sample is a control and `0` otherwise.

``` r
fit_data = fit_controls(analysis_path=analysis_path,
                        counts = filtered_counts,
                        sample_types=sample_types,
                        iterations = 20000,
                        burn_in = 5000,
                        precision_prior = c(1.5, 1000000) # shape and scale of gamma prior on the sample inverse dispersion parameter 
                        )
```

This returns a list including `map_means`, which is a vector of MAP estimates for the loci count proportion bias **m**. It also includes `run_data`, which is the loci counts for the samples we are interested in.

### Running the MCMC to infer the relative copy number profiles using `run_MCMC()`

Get data ready to run the MCMC.

``` r
# read cytoband information cytoBandIdeo.txt file - obtained from UCSC

cytoband_file = system.file("extdata",
                            "cytoBandIdeo.txt",
                            package = "conliga",
                            mustWork = TRUE)

cytoband = read_cytobands(cytoband_file,
                          chr_order=chromosomes)

cytoband
#> # A tibble: 811 x 5
#>    chr      start      end name   stain 
#>    <fct>    <int>    <int> <chr>  <chr> 
#>  1 chr1         0  2300000 p36.33 gneg  
#>  2 chr1   2300000  5300000 p36.32 gpos25
#>  3 chr1   5300000  7100000 p36.31 gneg  
#>  4 chr1   7100000  9100000 p36.23 gpos25
#>  5 chr1   9100000 12500000 p36.22 gneg  
#>  6 chr1  12500000 15900000 p36.21 gpos50
#>  7 chr1  15900000 20100000 p36.13 gneg  
#>  8 chr1  20100000 23600000 p36.12 gpos25
#>  9 chr1  23600000 27600000 p36.11 gneg  
#> 10 chr1  27600000 29900000 p35.3  gpos25
#> # ... with 801 more rows
# alternatively get cytoband information by querying UCSC

# cytoband = get_cytobands(genome="hg38",
#                          chr_order=chromosomes)

# we need to state our priors in the following format
priors_file = system.file("extdata",
                          "MCMC_priors.csv",
                          package = "conliga",
                          mustWork = TRUE)

priors = readr::read_csv(priors_file)

priors
#> # A tibble: 10 x 2
#>    prior                     value
#>    <chr>                     <dbl>
#>  1 hpp_rho_c              100000  
#>  2 hpp_rho_d                 100  
#>  3 hpp_g_a                     7  
#>  4 hpp_g_b                     1  
#>  5 hpp_apk_a                2000  
#>  6 hpp_apk_b                  10  
#>  7 precision_prior_shape       1.5
#>  8 precision_prior_scale 1000000  
#>  9 rcn_prior_shape             3  
#> 10 rcn_prior_scale             1

# which samples do we want to infer relative copy number profiles for?
samples_to_run = c("SLX-9396.FastSeqB",
                   "SLX-9396.FastSeqC")

# or run all non-controls:
# samples_to_run = sample_types %>%
#   dplyr::filter(Control == 0) %>%
#   .$Sample

num_cores = 2
iterations = 50000
thin = 10
```

Now we run the function `run_MCMC`.

``` r

run_MCMC(analysis_path=analysis_path,
         cytoband=cytoband,
         samples=samples_to_run,
         fit_data=fit_data,
         priors=priors,
         num_cores=num_cores,
         iterations=iterations,
         thin=thin,
         gamma = 1, # set gamma set to 1,
         max_states=30,
         init_states=TRUE) # start chain off by k-means clustering
```

Note that we have fixed the hyperparameter gamma to 1. The other hyperparameters alpha+kappa and rho will be inferred. The `run_MCMC` function usually takes approximately 30 mins per 10,000 iterations of the MCMC (depending on your machine).

### Processing the MCMC chains and obtaining results with `process_MCMC()`

Now we process the MCMC using the the `process_MCMC` function. Here we have selected to burn the first 5000 iterations.

``` r

burn_in = 5000

# samples %>% purrr::map(~process_MCMC(analysis_path=analysis_path,
#                                      sample_name=.x,
#                                      run_data=fit_data$run_data,
#                                      loci_means=fit_data$map_means,
#                                      cytoband=cytoband,
#                                      burn_in=burn_in))

samples %>% purrr::map(~process_MCMC(analysis_path=analysis_path,
                                     sample_name=.x,
                                     run_data=fit_data$run_data,
                                     loci_means=fit_data$map_means,
                                     cytoband=cytoband,
                                     burn_in=burn_in))
```

This function saves the inferred relative copy number calls to `analysis_path/results/sample_name/sample_name.tsv` and provides a number of basic plots under `analysis_path/results/sample_name/standard_plots/`.

Finally, we can delete the MCMC output files if we wish.

``` r
unlink(file.path(analysis_path, "MCMC"), recursive = TRUE, force = TRUE)
```
