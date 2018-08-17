sim_counts_from_LCN_with_different_precisions <- function(lociCN, chr_mat, mean_vec, total_draws,
                                                          precisions, trans_mat, init_dist, hidden_states) {
  return(mapply(simCountsFromLociCopyNum, precision=precisions,
                MoreArgs = list(lcns=lociCN,
                                chr_mat=chr_mat,
                                means=mean_vec,
                                total_draws=total_draws,
                                trans_mat=trans_mat,
                                init_dist=init_dist,
                                hs=hidden_states), SIMPLIFY=FALSE))
}

sim_counts_from_LCN_with_different_total_draws <- function(lociCN, chr_mat, mean_vec, total_draws_vec,
                                                           precision, trans_mat, init_dist, hidden_states) {
  return(mapply(simCountsFromLociCopyNum, total_draws=total_draws_vec,
                MoreArgs = list(lcns=lociCN,
                                chr_mat=chr_mat,
                                means=mean_vec,
                                precision=precisions,
                                trans_mat=trans_mat,
                                init_dist=init_dist,
                                hs=hidden_states), SIMPLIFY=FALSE))
}

sim_counts_from_LCN_with_different_mean_vecs <- function(lociCN, chr_mat, mean_vecs_list, total_draws,
                                                         precision, trans_mat, init_dist, hidden_states) {
  return(mapply(simCountsFromLociCopyNum, means=mean_vecs_list,
                MoreArgs = list(lcns=lociCN,
                                chr_mat=chr_mat,
                                precision=precisions,
                                total_draws=total_draws,
                                trans_mat=trans_mat,
                                init_dist=init_dist,
                                hs=hidden_states), SIMPLIFY=FALSE))
}
