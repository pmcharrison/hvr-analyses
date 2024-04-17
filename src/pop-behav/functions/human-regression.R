loadNamespace("R.utils")

run_human_regressions <- function(parent_dir, 
                                  out_dir,
                                  chosen_order, 
                                  max_eval = 1000L,
                                  optim_method = "NLOPT_LN_SBPLX",
                                  perm_int_reps = 5) {
  
  about <- readRDS(file.path(parent_dir, "about.rds"))
  
  run <- about$runs$ideal_regression %>% filter(ppm_order == chosen_order & 
                                                  viewpoint_group == "all")
  stopifnot(nrow(run) == 1)
  
  run_human_regression(run, parent_dir, out_dir, max_eval, about, optim_method, perm_int_reps)
}

run_human_regression <- function(run, 
                                 parent_dir,
                                 out_dir, 
                                 max_eval,
                                 about, 
                                 optim_method,
                                 perm_int_reps) {
  . <- get_data(parent_dir, run)
  par_scale <- pmax(1, abs(.$ideal_par))
  bounds <- get_bounds(.$predictors, .$ideal_par)
  ideal_model_cor <- get_ideal_model_cor(., par_scale, about)
  res <- optim_human(., par_scale, bounds, optim_method, max_eval, about)
  out <- make_human_regression_model(.,
                                     res,
                                     optim_method, 
                                     par_scale,
                                     bounds,
                                     ideal_model_cor,
                                     perm_int_reps)
  R.utils::mkdirs(out_dir)
  saveRDS(out, file.path(out_dir, "results.rds"))
}

get_data <- function(parent_dir, run) {
  stopifnot(nrow(run) == 1)
  
  viewpoint_dir <- file.path(parent_dir, "0-viewpoints")
  model_matrix_dir <- file.path(parent_dir, run$model_matrix_dir)
  ideal_regression_dir <- file.path(parent_dir, run$ideal_regression_dir)
  
  . <- list()
  
  .$ideal_mod <- readRDS(file.path(ideal_regression_dir, "results.rds"))
  .$ideal_par <- .$ideal_mod$par
  .$corpus <- readRDS(file.path(model_matrix_dir, "corpus.rds"))
  .$observation_matrix <- readRDS(file.path(model_matrix_dir, "observation-matrix.rds"))
  .$continuation_matrices <- readRDS(file.path(model_matrix_dir, "continuation-matrices.rds"))
  .$legal <- readRDS(file.path(model_matrix_dir, "legal.rds"))
  .$predictors <- readRDS(file.path(model_matrix_dir, "predictors.rds"))
  .$poly_degree <- readRDS(file.path(model_matrix_dir, "about.rds"))$poly_degree
  .$poly_coefs <-  readRDS(file.path(model_matrix_dir, "poly-coefs.rds"))
  .$moments <- readRDS(file.path(model_matrix_dir, "moments.rds"))
  .$viewpoint_labels <- readRDS(file.path(viewpoint_dir, "about.rds"))$viewpoint_labels
  
  return(.)
}

get_bounds <- function(predictors, ideal_par) {
  # higher-order polynomials are constrained to zero
  # discrete viewpoints are constrained to be non-negative
  stopifnot(all(names(ideal_par) == predictors$label))
  lower_bound <- if_else(predictors$discrete, 
                         0, 
                         if_else(predictors$poly_degree > 1, 0, -Inf)) 
  upper_bound <- if_else(!predictors$discrete & predictors$poly_degree > 1, 
                         0, Inf)
  list(lower = lower_bound,
       upper = upper_bound)
}

get_ideal_model_cor <- function(., par_scale, about) {
  res <- human_regression_eval_par(par = .$ideal_par / par_scale, 
                                   par_scale = par_scale, 
                                   observation_matrix = .$observation_matrix, 
                                   continuation_matrices = .$continuation_matrices, 
                                   legal = .$legal, 
                                   about = about, 
                                   corpus = .$corpus,
                                   counter = FALSE) * - 1
  message("Ideal observer correlation: ", round(res, digits = 5))
  reset_human_regression_eval_count()
  res
}

optim_human <- function(., par_scale, bounds, optim_method, max_eval, about) {
  withr::with_seed(1, {
    res <- nloptr::nloptr(
      x0 = rep(0, times = length(.$ideal_par)),
      eval_f = human_regression_eval_par,
      eval_grad_f = NULL,
      lb = bounds$lower,
      ub = bounds$upper,
      opts = list(algorithm = optim_method,
                  maxeval = max_eval,
                  ranseed = 1),
      par_scale = par_scale,
      about = about, 
      observation_matrix = .$observation_matrix, 
      continuation_matrices = .$continuation_matrices, 
      legal = .$legal,
      corpus = .$corpus,
      counter = TRUE
    )
    res$eval_f <- NULL 
    res$nloptr_environment <- NULL
    res
  })
}

make_human_regression_model <- function(.,
                                        res,
                                        optim_method, 
                                        par_scale,
                                        bounds,
                                        ideal_model_cor,
                                        perm_int_reps) {
  
  rescaled_par <- res$solution * par_scale
  
  mod <- hvr::new_regression_model(
    par = rescaled_par,
    cost = hvr:::cost(weights = rescaled_par, 
                      observation_matrix = .$observation_matrix,
                      continuation_matrices = .$continuation_matrices, 
                      legal = .$legal),
    corpus = .$corpus, 
    observation_matrix = .$observation_matrix, 
    continuation_matrices = .$continuation_matrices, 
    legal = .$legal, 
    predictors = .$predictors,
    perm_int = TRUE, 
    perm_int_seed = 1,
    perm_int_reps = perm_int_reps, 
    poly_degree = .$poly_degree, 
    poly_coefs = .$poly_coefs,
    moments = .$moments,
    viewpoint_labels = .$viewpoint_labels,
    optim_method = optim_method, 
    par_scale = par_scale,
    lower_bound = bounds$lower,
    optim_report = res
  )
  
  out <- list(mod = mod, 
              ideal_model_cor = ideal_model_cor, 
              human_model_cor = - res$objective,
              bounds = bounds)
  
  class(out) <- c("human_regression", class(out))
  
  out
}

human_regression_eval_count <- 0L
reset_human_regression_eval_count <- function() human_regression_eval_count <<- 0L

human_regression_eval_par <- function(par,
                                      par_scale,
                                      about,
                                      observation_matrix, 
                                      continuation_matrices, 
                                      legal,
                                      corpus,
                                      counter = TRUE) {
  df <- compare_regression_to_humans(par = par * par_scale,
                                     about, 
                                     observation_matrix, 
                                     continuation_matrices, 
                                     legal,
                                     corpus)
  c <- tryCatch(cor(df$information_content,
                    df$surprisal_rating_mean,
                    method = "spearman"),
                warning = function(...) 0)
  if (counter) human_regression_eval_count <<- human_regression_eval_count + 1L
  if (counter) message("i = ", human_regression_eval_count, ", rho = ", c)
  - c # the optimiser is a minimiser
}

compare_regression_to_humans <- function(
  par, 
  about = readRDS(file.path(parent_dir, "about.rds")), 
  observation_matrix = readRDS(file.path(model_matrix_dir, "observation-matrix.rds")), 
  continuation_matrices = readRDS(file.path(model_matrix_dir, "continuation-matrices.rds")), 
  legal = readRDS(file.path(model_matrix_dir, "legal.rds")),
  corpus = readRDS(file.path(model_matrix_dir, "corpus.rds")),
  parent_dir = stop("either 'about' or 'parent_dir' must be provided"),
  model_matrix_dir = stop("either model_matrix_dir must be provided, ",
                          "or observation_matrix, continuation_matrices, ",
                          "legal, and corpus must be provided individually")
) {
  corpus %>% 
    add_column(information_content = - log2(hvr:::event_probs(par, 
                                                             observation_matrix, 
                                                             continuation_matrices, 
                                                             legal))) %>% 
    rename(pos = event_id) %>% 
    mutate(stimulus_id = about$test_stimulus_ids[seq_id]) %>% 
    left_join(about$detail %>% select(stimulus_id, composition_id, event_id, pos,
                                      surprisal_rating_mean, 
                                      surprisal_rating_se,
                                      midi), 
              by = c("stimulus_id", "pos"))
}
