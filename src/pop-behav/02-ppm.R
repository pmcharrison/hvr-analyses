library(hrep)
library(tibble)

out_dir <- "output/pop-behav"
R.utils::mkdirs("1-ppm")

about <- readRDS(file.path(out_dir, "about.rds"))
runs <- about$runs$ppm

for (i in seq_len(nrow(runs))) {
  message("Performing PPM run ", i, " out of ", nrow(runs), "...")
  ppm_order <- runs$ppm_order[i]
  ppm_training_size <- runs$ppm_training_size[i]
  pretrain_ids <- runs$train_ids[[i]]
  run_dir <- file.path(out_dir, runs$ppm_dir[i])
  stopifnot(ppm_training_size == length(pretrain_ids))
  
  message(paste0("  training set size = ", ppm_training_size, 
                 ", order = ", ppm_order, 
                 ", output = ", run_dir))
  
  hvr::compute_ppm_analyses(
    parent_dir = out_dir,
    output_dir = run_dir,
    stm_opt = hvr::stm_options(order_bound = ppm_order),
    ltm_opt = hvr::ltm_options(order_bound = ppm_order),
    seq_pretrain = pretrain_ids
  )
}
