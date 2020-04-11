library(tidyverse)
theme_set(theme_classic() +
            ggplot2::theme(axis.text = ggplot2::element_text(colour = "black"),
                           axis.ticks = ggplot2::element_line(colour = "black")))

parent_dir <- "output/pop-behav"
about <- readRDS(file.path(parent_dir, "about.rds"))

out_dir <- file.path(parent_dir, "7-incon-replication")
R.utils::mkdirs(out_dir)

model_matrix_dir <- file.path(out_dir, "model-matrix")
regression_dir <- file.path(out_dir, "regression")

viewpoints <- c("har_18_harmonicity", "hutch_78_roughness", "num_pcs")

hvr::compute_model_matrix(
  parent_dir = parent_dir, 
  filter_corpus = function(x) x[x$event_id == 6, ], 
  poly_degree = 3, 
  ppm_dir = file.path(parent_dir, "runs/run-1/ppm"),
  output_dir = model_matrix_dir,
  viewpoints = viewpoints
)

mod <- hvr::viewpoint_regression(parent_dir = parent_dir, 
                                 model_matrix_dir = model_matrix_dir,
                                 output_dir = regression_dir,
                                 max_iter = 200,
                                 perm_int_reps = 5, 
                                 allow_negative_weights = FALSE)

viewpoint_labels <- mod$viewpoint_labels
hvr::plot_marginals(mod, viewpoint_labels = viewpoint_labels) + 
  theme(aspect.ratio = 1, panel.spacing = unit(7, "mm"))

ggsave(file.path(out_dir, "marginals.pdf"), width = 4.7, height = 4.7)

