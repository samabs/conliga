run_sticky <- function(counts,
                       chr_mat,
                       means,
                       max_states=30,
                       precision=500000,
                       iterations=10000,
                       jump_size_prob=c(0.8, 0.2),
                       gamma_a=3,
                       gamma_scale=1,
                       sample_ref="",
                       thin=1,
                       norm_sigma_s=0.15,
                       norm_sigma_l=0.5,
                       precision_sigma_s=50000,
                       precision_sigma_l=100000,
                       hpp_g_a=0.01,
                       hpp_g_b=0.01,
                       hpp_ak_a=0.01,
                       hpp_ak_b=0.01,
                       hpp_r_c=20,
                       hpp_r_d=1,
                       precision_p_shape=3.874839,
                       precision_p_scale = 84594.12,
                       gamma=-1,
                       rho=-1,
                       alpha_plus_kappa=-1,
                       init_hidden_states,
                       init_state_rcns) {

  if (missing(init_hidden_states)) {
    init_hidden_states = rep(0, length(counts))
  }
  
  if (missing(init_state_rcns)) {
    init_state_rcns = rep(0, max_states)
  }
  
  runStickyHMM(counts,
               chr_mat,
               max_states,
               means,
               precision,
               iterations,
               jump_size_prob,
               gamma_a,
               gamma_scale,
               sample_ref,
               thin,
               norm_sigma_s,
               norm_sigma_l,
               precision_sigma_s,
               precision_sigma_l,
               hpp_g_a,
               hpp_g_b,
               hpp_ak_a,
               hpp_ak_b,
               hpp_r_c,
               hpp_r_d,
               precision_p_shape,
               precision_p_scale,
               gg=gamma,
               rr=rho,
               apk=alpha_plus_kappa,
               init_hidden_states=init_hidden_states,
               init_state_rcn=init_state_rcns)
}

run_sticky_in_parallel <- function(num_cores,
                                   counts_list,
                                   chr_mat_list,
                                   means_list, 
                                   hpp_g_a_list,
                                   hpp_g_b_list,
                                   hpp_ak_a_list,
                                   hpp_ak_b_list, 
                                   hpp_r_c_list,
                                   hpp_r_d_list,
                                   sample_ref_list, 
                                   precision_p_shape_list,
                                   precision_p_scale_list,
                                   rcn_shape_list,
                                   rcn_scale_list,
                                   gamma_list,
                                   rho_list,
                                   apk_list,
                                   init_hidden_states_list,
                                   init_state_rcns_list,
                                   max_states=30,
                                   precision=500000,
                                   iterations=10000,
                                   jump_size_prob=c(0.8, 0.2),
                                   thin=1,
                                   norm_sigma_s=0.15,
                                   norm_sigma_l=0.5,
                                   precision_sigma_s=50000,
                                   precision_sigma_l=100000) {
  
  if (missing(init_hidden_states_list)) {
    rep(list(rep(0, length(c(1,2,3,4)))), 24)
    init_hidden_states_list = rep(list(rep(0, length(counts_list[[1]]))), length(counts_list))
  }
  
  if (missing(init_state_rcns_list)) {
    init_state_rcns_list = rep(list(rep(0, max_states)), length(counts_list))
  }
  
  return(system.time(
    parallel::mcmapply(FUN=runStickyHMM,
                       counts=counts_list,
                       sample_ref=sample_ref_list,
                       chr_mat=chr_mat_list,
                       means=means_list,
                       hpp_g_a=hpp_g_a_list,
                       hpp_g_b=hpp_g_b_list,
                       hpp_ak_a=hpp_ak_a_list,
                       hpp_ak_b=hpp_ak_b_list,
                       hpp_r_c=hpp_r_c_list,
                       hpp_r_d=hpp_r_d_list,
                       precision_p_shape=precision_p_shape_list,
                       precision_p_scale=precision_p_scale_list,
                       gamma_a=rcn_shape_list,
                       gamma_scale=rcn_scale_list,
                       gg=gamma_list,
                       rr=rho_list,
                       apk=apk_list,
                       init_hidden_states=init_hidden_states_list,
                       init_state_rcn=init_state_rcns_list,
                       mc.preschedule = FALSE,
                       mc.cores = num_cores,
                       MoreArgs=list(
                         norm_sigma_s=norm_sigma_s,
                         norm_sigma_l=norm_sigma_l,
                         precision_sigma_s=precision_sigma_s,
                         precision_sigma_l=precision_sigma_l,
                         max_states=max_states,
                         precision=precision,
                         iterations=iterations,
                         jump_size_prob=jump_size_prob,
                         thin=thin)
    )
  ))
}

