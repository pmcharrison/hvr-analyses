library(hrep)

out_dir <- "output/pop-behav"
about <- readRDS(file.path(out_dir, "about.rds"))

hvr::compute_viewpoints(
  corpus = about$combined_corpus,
  parent_dir = out_dir,
  seq_test = about$indices$test,
  viewpoints = hvr::hvr_viewpoints
)
