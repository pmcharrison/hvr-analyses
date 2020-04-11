# library(tidyverse)
# theme_set(theme_classic() +
#             ggplot2::theme(axis.text = ggplot2::element_text(colour = "black"),
#                            axis.ticks = ggplot2::element_line(colour = "black")))
# source("src/pop-behav/functions/human-regression.R")
# 
# target_chord <- 6L
# 
# parent_dir <- "output/pop-behav"
# about <- readRDS(file.path(parent_dir, "about.rds"))
# viewpoint_dir <- file.path(parent_dir, "0-viewpoints")
# 
# ppm_order <- 1
# ppm_training_size <- 403 
# 
# run_dir <- about$runs %>% filter(ppm_order == !!ppm_order,
#                                  ppm_training_size == !!ppm_training_size) %>%
#   pull(dir) %>%
#   file.path(parent_dir, "runs", .)
# 
# spec <- list(
#   # continuous_only = c("har_18_harmonicity", "hutch_78_roughness", "spec_dist_3", "vl_dist", "num_pcs"),
#   # just_vl_dist = "vl_dist",
#   # v3 = c("vl_dist", "num_pcs", "hutch_78_roughness")
#   # continuous_only_pcs = c("har_18_harmonicity", "hutch_78_roughness", "spec_dist_3", "vl_dist", "num_pcs")
#   discrete_only = readRDS(file.path(viewpoint_dir, "about.rds"))$discrete_viewpoints %>%
#     map_chr("name") %>% unname()
# )
# 
# map2(spec, names(spec), function(viewpoints, run_label) {
#   message("Conducting run: ", run_label)
#   message("  viewpoints: ", paste(viewpoints, collapse = ", "))
#   
#   out_dir <- file.path(parent_dir, "4-benchmark-models", run_label)
#   R.utils::mkdirs(out_dir)
#   
#   model_matrix_dir <- file.path(out_dir, "model-matrix")
#   regression_dir <- file.path(out_dir, "regression")
#   
#   message("  output directory: ", out_dir)
# 
#   hvr::compute_model_matrix(
#     parent_dir = parent_dir, 
#     filter_corpus = function(x) x[x$event_id == 6, ], 
#     poly_degree = 4, 
#     ppm_dir = file.path(run_dir, "ppm"),
#     output_dir = model_matrix_dir,
#     viewpoints = viewpoints
#   )
#   
#   mod <- hvr::viewpoint_regression(parent_dir = parent_dir, 
#                                    model_matrix_dir = model_matrix_dir,
#                                    output_dir = regression_dir,
#                                    max_iter = 100,
#                                    perm_int_reps = 5, 
#                                    allow_negative_weights = FALSE)
#   
#   compare <- compare_regression_to_humans(mod$par, parent_dir = parent_dir, model_matrix_dir = model_matrix_dir)
#   write_csv(compare, file.path(out_dir, "human-cor.csv"))
#   
#   "Spearman correlation between model output and mean surprisal ratings: %s" %>% 
#     sprintf(cor(compare$surprisal_rating_mean, compare$information_content, method = "spearman")) %>% 
#     writeLines(file.path(out_dir, "human-cor.txt"))
#   
#   ggplot(compare, aes(information_content, surprisal_rating_mean)) +
#     geom_point(shape = 21) +
#     scale_x_continuous("Information content") +
#     scale_y_continuous("Mean surprisal rating") +
#     theme(aspect.ratio = 1)
#   ggsave(file.path(out_dir, "human-cor.pdf"), width = 3.5, height = 3.5)
#   
#   # We could add a human regression model, but maybe it's not worth the effort.
#   # If we do so, we'd need to refactor the human regression section.
#   
#   viewpoint_labels <- mod$viewpoint_labels %>% 
#     mutate(viewpoint_label = recode(viewpoint_label, 
#                                     `Spectral distance (HL = 3)` = "Spectral distance"))
#   hvr::plot_marginals(mod, viewpoint_labels = viewpoint_labels, scales = "free") + 
#     theme(aspect.ratio = 1, panel.spacing = unit(7, "mm"))
#   
#   ggsave(file.path(out_dir, "marginals.pdf"), width = 4.7, height = 4.7)
#   message("Completed run.\n")
# })
