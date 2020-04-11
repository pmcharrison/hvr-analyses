plot_discrete_viewpoint <- function(viewpoint_sequences, 
                                    viewpoint,
                                    alphabet_size,
                                    zero_indexed = FALSE) {
  checkmate::qassert(viewpoint, "S1")
  checkmate::qassert(alphabet_size, "X1")
  
  seqs <- viewpoint_sequences %>% 
    purrr::map("discrete") %>% 
    purrr::map(~ .[viewpoint, ]) %>% 
    purrr::map(na.omit)
  
  mod <- ppm::new_ppm_simple(alphabet_size = alphabet_size, 
                             order_bound = 1)
  for (s in seqs) ppm::model_seq(mod, s, predict = FALSE)
  
  ppm::plot_n_grams(mod,
                    zero_indexed = zero_indexed,
                    bigram_fill_scale = ggplot2::scale_fill_viridis_c("Probability (relative)"))
}
