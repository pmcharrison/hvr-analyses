check_df <- function(df, n_sample = 100) {
  stopifnot(n_sample <= nrow(df))
  rows <- sample(nrow(df), n_sample)
  for (i in seq_along(rows)) {
    check_row(df[rows[i], ])
  }
  df
}

only_alpha <- function(x) {
  gsub("[^a-zA-Z]", "", x)
}

check_row <- function(row) {
  # Check composition ID versus title
  hcorp::popular_1[[row$composition_id]] %>% 
    metadata() %>% 
    `$`(description) %>% 
    only_alpha() %>% 
    expect_equal(row$title %>% only_alpha())
  
  # Check chord
  hcorp::popular_1[[row$composition_id]] %>% 
    decode() %>% 
    `[[`(row$event_id) %>% 
    expect_equal(pc_chord(row$midi))
}

get_stimulus_data <- function(file) {
  read_csv(file, col_types = cols()) %>% 
    rename(pos = chord_id) %>% 
    select(- dataset_id) %>% 
    mutate(composition_id = as.integer(composition_id) + 1L,
           event_id = as.integer(event_id) + 1L) %>% 
    arrange(stimulus_id, pos) %>% 
    check_df()
}

get_test_corpus <- function(stimulus_ids, df) {
  stimulus_ids %>% 
    map(function(id) {
      x <- filter(df, stimulus_id == id)
      composition_id <- unique(x$composition_id)
      expect_equal(x$pos, 1:8)
      stopifnot(length(composition_id) == 1)
      hcorp::popular_1[[composition_id]][x$event_id]
    }) %>% 
    hrep::corpus(type = "pc_chord")
}

get_train_corpus <- function(df) {
  train_compositions <- setdiff(seq_along(hcorp::popular_1), 
                                df$composition_id %>% unique() %>% sort())
  x <- hcorp::popular_1[train_compositions]
  metadata(x) <- list()
  x
}

save_about <- function(out_dir, df, test_corpus, train_corpus, combined_corpus, 
                       test_ids, train_ids, stimulus_ids, runs) {
  about <- list(
    detail = df,
    test_corpus = test_corpus, # for reference only
    train_corpus = train_corpus, # for reference only
    combined_corpus = combined_corpus, # this goes into the hvr modelling
    indices = list(test = test_ids, # which parts of combined_corpus are test or training?
                   train = train_ids),
    test_stimulus_ids = stimulus_ids, # in the order of indices$test
    runs = runs
  )
  R.utils::mkdirs(out_dir)
  saveRDS(about, file.path(out_dir, "about.rds"))
}

get_runs <- function(ppm_orders, ppm_training_sizes, poly_degrees, all_train_ids, viewpoint_group) {
  ppm_runs <- expand.grid(ppm_order = ppm_orders, 
                          ppm_training_size = ppm_training_sizes) %>% 
    as_tibble() %>% 
    mutate(ppm_run_id = seq_along(ppm_order),
           ppm_dir = paste("1-ppm/run-", ppm_run_id, sep = "")) %>% 
    left_join(., sample_train_ids(.$ppm_training_size, all_train_ids),
              by = "ppm_training_size") %>% 
    arrange(ppm_run_id) %>% 
    select(ppm_run_id, ppm_dir, everything()) %>% 
    check_ppm_runs(ppm_orders, ppm_training_sizes, all_train_ids)
  
  model_matrix_runs <- expand.grid(ppm_run_id = ppm_runs$ppm_run_id,
                                   viewpoint_group = viewpoint_group,
                                   ltm = c(TRUE, FALSE),
                                   poly_degree = poly_degrees,
                                   stringsAsFactors = FALSE) %>% 
    rbind(expand.grid(ppm_run_id = ppm_runs$ppm_run_id,
                      viewpoint_group = "pc_chord_only",
                      ltm = c(TRUE, FALSE),
                      poly_degree = poly_degrees,
                      stringsAsFactors = FALSE)) %>% 
    as_tibble() %>% 
    left_join(ppm_runs, by = "ppm_run_id") %>% 
    filter(viewpoint_group != "continuous_only" | ppm_order == 0) %>% 
    filter(!(viewpoint_group == "continuous_only" & ltm)) %>% 
    mutate(model_matrix_run_id = seq_along(viewpoint_group),
           model_matrix_dir = paste("2-model-matrix/run-", model_matrix_run_id, sep = "")) %>% 
    select(model_matrix_run_id, model_matrix_dir, everything())
  
  ideal_regression_runs <- model_matrix_runs %>% 
    mutate(ideal_regression_run_id = seq_along(model_matrix_run_id),
           ideal_regression_dir = paste("3-ideal-regression/run-",
                                        ideal_regression_run_id, 
                                        sep = "")) %>% 
    select(ideal_regression_run_id, ideal_regression_dir, everything())
  
  list(ppm = ppm_runs,
       model_matrix = model_matrix_runs,
       ideal_regression = ideal_regression_runs)
}

check_ppm_runs <- function(runs, ppm_orders, ppm_training_sizes, all_train_ids) {
  stopifnot(all(runs$ppm_run_id == seq_len(nrow(runs))),
            !anyDuplicated(select(runs, ppm_order, ppm_training_size)))
  sizes <- sort(unique(runs$ppm_training_size))
  
  # Check that each training set includes the smaller training sets,
  # and that for a given training set size all training sets are identical.
  for (i in seq_along(sizes)) {
    if (i > 1) {
      ids <- map(c(sizes[i - 1], sizes[i]), function(size) {
        x <- runs %>% filter(ppm_training_size == size) %>% pull(train_ids) %>% unique()
        stopifnot(length(x) == 1)
        x[[1]]
      })
      stopifnot(all(ids[[1]] %in% ids[[2]]))
    }
  }
  runs
}

sample_train_ids <- function(ppm_training_sizes, all_train_ids) {
  withr::with_seed(1, {
    unique_sizes <- sort(unique(ppm_training_sizes))
    max_sample <- all_train_ids[sample(x = length(all_train_ids),
                                       size = max(unique_sizes),
                                       replace = FALSE)]
    tibble(
      ppm_training_size = unique_sizes,
      train_ids =  map(unique_sizes, ~ sort(max_sample[seq_len(.)]))
    )
  })
}
