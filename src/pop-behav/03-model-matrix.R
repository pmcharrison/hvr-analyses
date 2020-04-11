library(hrep)
library(tibble)

out_dir <- "output/pop-behav"
R.utils::mkdirs(file.path(out_dir, "2-model-matrix"))

about <- readRDS(file.path(out_dir, "about.rds"))
runs <- about$runs$model_matrix

viewpoint_groups <- list(
  all = hvr::hvr_viewpoints %>%
    names() %>% 
    setdiff("num_pcs"),
  continuous_only = hvr::hvr_viewpoints %>%
    Filter(Negate(hvr::is_discrete), .) %>% 
    names() %>% 
    setdiff("num_pcs"),
  categorical_only = hvr::hvr_viewpoints %>% 
    Filter(hvr::is_discrete, .) %>% 
    names(),
  pc_chord_only = "pc_chord"
)

for (i in seq_len(nrow(runs))) {
  message("Creating model matrix ", i, " out of ", nrow(runs), "...")
  ppm_dir <- file.path(out_dir, runs$ppm_dir[i])
  model_matrix_dir <- file.path(out_dir, runs$model_matrix_dir[i])
  
  viewpoint_group <- runs$viewpoint_group[i]
  viewpoints <- viewpoint_groups[[viewpoint_group]]
  stopifnot(is.character(viewpoints))
  
  ltm <- runs$ltm[i]
  
  message("Analysing viewpoint group '", viewpoint_group, "': ",
          paste(viewpoints, collapse = ", "), ".")
  message("LTM = ", ltm)
  
  R.utils::mkdirs(model_matrix_dir)
  
  hvr::compute_model_matrix(
    parent_dir = out_dir, 
    filter_corpus = function(x) x[x$event_id == 6, ], 
    poly_degree = runs$poly_degree[i], 
    ppm_dir = ppm_dir, 
    output_dir = model_matrix_dir,
    viewpoints = viewpoints %>% setNames(., .), 
    ltm = ltm
  )
} 
