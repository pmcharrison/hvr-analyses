# This script takes about 12 hours to complete 
# (the latter runs are faster than the earlier runs)

library(hrep)
library(tibble)

out_dir <- "output/pop-behav"
R.utils::mkdirs(file.path(out_dir, "3-ideal-regression"))

about <- readRDS(file.path(out_dir, "about.rds"))
runs <- about$runs$ideal_regression

for (i in seq_len(nrow(runs))) {
  message("Conducting ideal viewpoint regression ", i, " out of ", nrow(runs), "...")
  model_matrix_dir <- file.path(out_dir, runs$model_matrix_dir[i])
  ideal_regression_dir <- file.path(out_dir, runs$ideal_regression_dir[i])
  message("  input directory = ", model_matrix_dir)
  message("  output directory = ", ideal_regression_dir, "\n")
  hvr::viewpoint_regression(parent_dir = out_dir, 
                            model_matrix_dir = model_matrix_dir,
                            output_dir = ideal_regression_dir,
                            max_iter = 500,
                            perm_int_reps = 1, 
                            allow_negative_weights = FALSE)
  gc()
} 

# viewpoint_labels <- tibble(viewpoint = purrr::map_chr(hvr::hvr_viewpoints, "name"),
#                            viewpoint_label = purrr::map_chr(hvr::hvr_viewpoints, "label"))
# res <- readRDS("output/pop-behav/runs/run_4/ideal-regression/results.rds")
# res$viewpoint_labels <- viewpoint_labels
# saveRDS(res, "output/pop-behav/runs/run_4/ideal-regression/results.rds")
# 
# hvr::plot_marginals(res)
