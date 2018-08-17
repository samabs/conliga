swap_labels <- function(MCMC_data, perm, burned_n_states, nn, n_states, iter_w_data, burn_in, map_num_states=TRUE) {
  iterations = MCMC_data$params$iterations
  thin = MCMC_data$params$thin
  if (map_num_states) {
    map_state = Mode(burned_n_states)
    uu = data.frame(MCMC_data$stateRCN, iter = iter_w_data)
    colnames(uu) = c(0:(ncol(MCMC_data$stateRCN)-1), "iter")
    uu = uu %>% dplyr::filter(iter > burn_in)
    
    uu = uu %>% dplyr::left_join(nn) %>%
      dplyr::filter(num_states == map_state)
    uu = uu %>% dplyr::select(-iter, -num_states)
    uu = uu %>% as.matrix()
    state_labels = matrix(NA, nrow=nrow(uu), ncol=map_state)
    subset_uu <- matrix(NA, nrow=nrow(uu), ncol=map_state)
    iter_hs = MCMC_data$hs %>%
      dplyr::mutate(iter = iter_w_data) %>%
      dplyr::filter(iter > burn_in)
    iter_hs = iter_hs %>%
      dplyr::left_join(nn) %>%
      dplyr::filter(num_states == map_state)
    iter_hs = iter_hs %>%
      dplyr::select(-iter, -num_states)
    
    iters = nn %>% dplyr::filter(num_states == map_state) %>% .[["iter"]]
    used_states = map_state
    
    
  } else {
    uu = data.frame(MCMC_data$stateRCN, iter = iter_w_data)
    colnames(uu) = c(0:(ncol(MCMC_data$stateRCN)-1), "iter")
    uu = uu %>%
      dplyr::filter(iter > burn_in)
    uu = uu %>% 
      dplyr::select(-iter)
    uu = uu %>% as.matrix()
    max_states = max(burned_n_states)
    state_labels = matrix(NA, nrow=nrow(uu), ncol=max_states)
    subset_uu <- matrix(NA, nrow=nrow(uu), ncol=max_states)
    iter_hs = MCMC_data$hs %>%
      dplyr::mutate(iter = iter_w_data) %>%
      dplyr::filter(iter > burn_in)
    iter_hs = iter_hs %>%
      dplyr::select(-iter)
    iters = nn %>% .[["iter"]]
    used_states = max_states
  }
  
  for(i in 1:nrow(uu)) {
    states <- unique(as.integer(iter_hs[i,]))
    states <- sort(states)
    state_labels[i,1:length(states)] <- states
    subset_uu[i,] <- uu[i,state_labels[i,]+1]
  }
  
  new_uu = subset_uu
  
  for (i in 1:nrow(perm)) {
    new_uu[i,] = subset_uu[i,perm[i,]+1]
  }
  
  dd = state_labels
  new_dd = dd
  
  for (i in 1:nrow(perm)) {
    new_dd[i,] = dd[i,perm[i,]+1]
  }

  uu = data.frame(new_uu, iter=iters)
  colnames(uu) = c(0:(used_states-1), "iter")
  uu = uu %>%
    tidyr::gather(state, CN, -iter)
  
  dd = data.frame(new_dd, iter=iters)
  colnames(dd) = c(0:(used_states-1), "iter")
  dd = dd %>%
    tidyr::gather(state, index, -iter)

  
  tt = uu %>% dplyr::left_join(dd) %>%
    dplyr::mutate(index = as.character(index))

  return(tt)
  
}