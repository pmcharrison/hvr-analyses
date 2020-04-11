get_ideal_model_predictions <- function(parent_dir, about, max_information_content) {
  message("Generating predictions from ideal-observer models...")
  about$runs$ideal_regression %>% 
    mutate(
      mod = map(ideal_regression_dir,
                ~  readRDS(file.path(parent_dir, ., "results.rds"))),
      predictions = plyr::llply(seq_along(mod), function(i) {
        compare_regression_to_humans(
          mod[[i]]$par,
          about = about,
          model_matrix_dir = file.path(parent_dir, model_matrix_dir[i])
        ) %>%
          mutate(information_content = pmin(information_content,
                                            max_information_content))
      }, .progress = "time"),
      spearman = map_dbl(predictions, 
                         ~ suppressWarnings(cor(.$information_content, 
                                                .$surprisal_rating_mean,
                                                method = "spearman"))),
      pearson = map(predictions,
                    ~ cor.test(.$information_content, 
                               .$surprisal_rating_mean,
                               method = "pearson") %>% 
                      broom::glance()),
      pearson_cor_att = map_dbl(predictions, 
                                ~ cor_att(.$information_content,
                                          .$surprisal_rating_mean,
                                          x_se = 0,
                                          y_se = mean(.$surprisal_rating_se))),
      adj_r2 = map_dbl(predictions, 
                       ~ get_adj_r2(.$information_content,
                                    .$surprisal_rating_mean)),
      n_cases = map_int(predictions, nrow),
    ) %>% 
    mutate(
      pearson_lower = map_dbl(pearson, "conf.low"),
      pearson_upper = map_dbl(pearson, "conf.high"),
      pearson = map_dbl(pearson, "estimate")
    ) %>% 
    select(- train_ids) %>% 
    add_cross_entropy()
}

save_ideal_model_predictions <- function(x, out_dir) {
  saveRDS(x, file.path(out_dir, "ideal-model-predictions.rds"))
  
  x %>% 
    select(ideal_regression_run_id, viewpoint_group, ltm, poly_degree, ppm_order,
           ppm_training_size, spearman,
           starts_with("pearson"), starts_with("cross_entropy")) %>% 
    write_csv(file.path(out_dir, "ideal-model-predictions.csv"))
  
  x
}

regress_models <- function(i, j, ideal_model_predictions) {
  get_dual_predictions(i, j, ideal_model_predictions) %>% 
    lm(scale(surprisal_rating_mean) ~ 
         scale(information_content.x) + 
         scale(information_content.y), 
       data = .)
}

get_benchmarks <- function(file = "input/pop-behav/pop-behav.csv") {
  message("Getting benchmarks...")
  x <- read_csv(file, col_types = cols()) %>%
    mutate(composition_id = composition_id + 1L,
           event_id = event_id + 1L) %>% 
    group_by(stimulus_id) %>% 
    mutate(full_midi = lapply(midi, function(...) midi)) %>%
    ungroup() %>% 
    na.omit() %>% 
    mutate(., 
           new_benchmarks = plyr::mlply(., get_new_benchmarks, .progress = "time")) %>% 
    select(composition_id, event_id, 
           surprisal_rating_mean,
           surprisal_rating_se,
           starts_with("leman"), "collins_2014", "incon",
           new_benchmarks)
  
  y <- bind_cols(x %>% select(- new_benchmarks),
                 x %>% pull(new_benchmarks) %>% bind_rows())
  
  attr(y, "benchmarks") <- setdiff(names(y), 
                                   c("composition_id", "event_id", 
                                     "surprisal_rating_mean",
                                     "surprisal_rating_se"))
  
  y
}

get_adj_r2 <- function(x, y) {
  tibble(x, y) %>% 
    lm(x ~ y, data = .) %>%
    broom::glance() %>% 
    pull(adj.r.squared)
}

eval_benchmarks <- function(benchmarks) {
  chosen_benchmarks <- c("leman_2000_v6", "collins_2014", "incon",
                         "parn94_pitch_commonality", "parn94_pitch_distance",
                         "pi_dist", "pc_dist", "spec_sim_4")
  
  x <- tibble(
    benchmark = attr(benchmarks, "benchmarks"),
    chosen = benchmark %in% chosen_benchmarks,
    spearman = map_dbl(benchmark, ~ cor(benchmarks$surprisal_rating_mean, 
                                        benchmarks[[.]],
                                        method = "spearman")),
    pearson = map(benchmark, ~ cor.test(benchmarks$surprisal_rating_mean, 
                                        benchmarks[[.]]) %>% broom::glance()),
    pearson_cor_att = map_dbl(benchmark, 
                              ~ cor_att(benchmarks$surprisal_rating_mean,
                                        benchmarks[[.]],
                                        mean(benchmarks$surprisal_rating_se),
                                        0)),
    adj_r2 = map_dbl(benchmark, ~ get_adj_r2(benchmarks$surprisal_rating_mean, 
                                             benchmarks[[.]]))
  ) %>% 
    mutate(
      pearson_lower = map_dbl(pearson, "conf.low"),
      pearson_upper = map_dbl(pearson, "conf.high"),
      pearson = map_dbl(pearson, "estimate")
    )
  
  combined_benchmark_mod <- chosen_benchmarks %>% 
    paste(collapse = " + ") %>% 
    paste0("surprisal_rating_mean ~ ", .) %>% 
    lm(data = benchmarks)
  
  combined_benchmark_data <- tibble(model = predict(combined_benchmark_mod, 
                                                    newdata = benchmarks),
                                    listener = benchmarks$surprisal_rating_mean,
                                    listener_se = benchmarks$surprisal_rating_se)
  
  y <- tibble(benchmark = "combined_benchmarks",
              spearman = cor(combined_benchmark_data$model, 
                             combined_benchmark_data$listener,
                             method = "spearman"),
              pearson = list(cor.test(combined_benchmark_data$model, 
                                      combined_benchmark_data$listener,
                                      method = "pearson") %>% broom::glance()),
              pearson_cor_att = cor_att(combined_benchmark_data$model,
                                        combined_benchmark_data$listener,
                                        0,
                                        mean(combined_benchmark_data$listener_se)),
              adj_r2 = combined_benchmark_mod %>% broom::glance() %>% pull(adj.r.squared)) %>% 
    mutate(
      pearson_lower = map_dbl(pearson, "conf.low"),
      pearson_upper = map_dbl(pearson, "conf.high"),
      pearson = map_dbl(pearson, "estimate")
    )
  
  z <- bind_rows(x, y)
  attr(z, "chosen_benchmarks") <- chosen_benchmarks
  attr(z, "combined_benchmark_data") <- combined_benchmark_data
  z
}

get_eval_ideal_vs_benchmarks <- function(benchmarks_eval, ideal_model_predictions) {
  x <- ideal_model_predictions %>% 
    filter(ppm_order == 5 & 
             viewpoint_group == "all" & 
             ltm) %>% 
    select(spearman, pearson, pearson_cor_att, adj_r2, 
           pearson_lower, pearson_upper) %>% 
    add_column(model = "ideal_observer", .before = 1)
  
  y <- benchmarks_eval %>% 
    filter(chosen | benchmark == "combined_benchmarks") %>% 
    select(- chosen) %>% 
    rename(model = benchmark)
  
  bind_rows(x, y)
}


get_new_benchmarks <- function(composition_id, event_id, chord_id, midi, full_midi, ...) {
  stim <- hcorp::popular_1[[composition_id]][seq(to = event_id, length.out = chord_id)] %>% 
    hrep::decode()
  
  full_midi <- full_midi[[1]][1:chord_id] %>% map(hrep::pi_chord)
  
  stopifnot(identical(hrep::pc_chord(midi),
                      hrep::pc_chord(stim[[chord_id]])),
            identical(map(full_midi, hrep::pc_chord),
                      stim[1:chord_id] %>% as.list()))
  
  final_pi_chords <- full_midi[seq(to = chord_id, 
                                   length = 2)] %>% set_names("x", "y")
  
  final_pitch_saliences <- map(final_pi_chords, parn94::pitch_salience)
  
  list(
    parn94_pitch_commonality = do.call(parn94::pitch_commonality, final_pitch_saliences),
    parn94_pitch_distance = do.call(parn94::pitch_distance, final_pitch_saliences),
    pi_dist = minVL::min_vl_dist(final_pi_chords[[1]], 
                                 final_pi_chords[[2]],
                                 elt_type = "pitch"),
    pc_dist = minVL::min_vl_dist(final_pi_chords[[1]] %>% hrep::pc_set(), 
                                 final_pi_chords[[2]] %>% hrep::pc_set(),
                                 elt_type = "pc"),
    spec_sim_ = 1:4 %>% map_dbl(~ specdec::spectral_similarity(stim, half_life = .)[chord_id])
    
  ) %>% 
    unlist() %>% 
    as.list() %>% 
    as_tibble()
}

benchmark_model <- function(i, ideal_model_predictions, file) {
  behav_data <- get_benchmarks()
  
  our_model <- ideal_model_predictions$predictions[[i]] %>% 
    select(composition_id, event_id, information_content)
  
  stopifnot(all.equal(behav_data$composition_id, 
                      our_model$composition_id),
            all.equal(behav_data$event_id,
                      our_model$event_id))
  
  all_data <- inner_join(behav_data, 
                         our_model, 
                         by = c("composition_id", "event_id"))
  
  . <- list()
  
  .$cors <- all_data %>% 
    select("surprisal_rating_mean", 
           "information_content", 
           starts_with("leman"), "collins_2014", "incon") %>% 
    cor() %>% 
    as_tibble() %>% 
    slice(1) %>% 
    select(- surprisal_rating_mean)
  
  .$lm_1 <- lm(surprisal_rating_mean ~ leman_2000_v6 + collins_2014 + incon,
               data = all_data)
  
  .$lm_2 <- update(.$lm_1, ~ . + information_content)
  
  .$anova <- anova(.$lm_1, .$lm_2)
  
  .
}

get_predictions <- function(i, ideal_model_predictions) {
  ideal_model_predictions$predictions[[i]] %>% 
    select(seq_event_id, information_content, surprisal_rating_mean)
}

get_dual_predictions <- function(i, j, ideal_model_predictions) {
  c(i, j) %>%
    map(get_predictions, ideal_model_predictions) %>% 
    {inner_join(.[[1]], .[[2]] %>% select(- surprisal_rating_mean),
                by = "seq_event_id")}
}

cor_model <- function(i, ideal_model_predictions) {
  ideal_model_predictions %>% 
    mutate(pearson_r_lower = map_dbl(pearson_cor, "conf.low"),
           pearson_r_upper = map_dbl(pearson_cor, "conf.high")) %>% 
    select(starts_with("pearson_r"), spearman_rho) %>% 
    slice(i) %>% 
    mutate_all(round, 2)
}

cor_models <- function(i, j, ideal_model_predictions, method = "pearson") {
  get_dual_predictions(i, j, ideal_model_predictions) %>% 
    with(cor.test(information_content.x, information_content.y, method = method))
}

plot_model <- function(i, ideal_model_predictions) {
  ideal_model_predictions$predictions[[i]] %>% 
    with(plot(information_content, surprisal_rating_mean))
}

get_cor_compare <- function(ideal_model_predictions) {
  n <- nrow(ideal_model_predictions)
  expand.grid(i = seq_len(n), 
              j = seq_len(n)) %>% 
    as_tibble() %>% 
    mutate(r_i = map_dbl(i, ~ ideal_model_predictions$pearson[.]),
           r_j = map_dbl(j, ~ ideal_model_predictions$pearson[.]),
           diff = r_i - r_j,
           conf = map2(i, j, compare_cors, ideal_model_predictions),
           conf_lower = map_dbl(conf, 1),
           conf_upper = map_dbl(conf, 2),
           sig = !(conf_lower <= 0 & conf_upper >= 0)) %>% 
    select(- conf)
}

compare_cors <- function(k, h, ideal_model_predictions) {
  data <- map(c(k, h), 
              ~ ideal_model_predictions$predictions[[.]] %>% 
                select(seq_event_id, information_content, surprisal_rating_mean)) %>% 
    {inner_join(.[[1]] %>% select(- surprisal_rating_mean), 
                .[[2]], 
                by = "seq_event_id", suffix = c("_k", "_h"))}
  
  r_jk <- cor(data$surprisal_rating_mean, 
              data$information_content_k)
  
  r_jh <- cor(data$surprisal_rating_mean, 
              data$information_content_h)
  
  r_kh <- cor(data$information_content_k, 
              data$information_content_h)
  
  test <- cocor::cocor.dep.groups.overlap(r.jk = r_jk, r.jh = r_jh, r.kh = r_kh,
                                          n = nrow(data), test = "zou2007")
  
  stopifnot(test@r.jk == ideal_model_predictions$pearson_r[k],
            test@r.jh == ideal_model_predictions$pearson_r[h])
  
  test@zou2007$conf.int
}

add_cross_entropy <- function(x) {
  x %>% 
    mutate(cross_entropy = map(predictions, ~ {
      tibble(
        mean = mean(.$information_content),
        sd = sd(.$information_content),
        se = sd / sqrt(nrow(.)),
        lower = mean - 1.96 * se,
        upper = mean + 1.96 * se
      ) %>% as.list()
    })) %>% 
    mutate(cross_entropy_lower = map_dbl(cross_entropy, "lower"),
           cross_entropy_upper = map_dbl(cross_entropy, "upper"),
           cross_entropy = map_dbl(cross_entropy, "mean"))
}

drop_l0 <- function(x) {
  gsub("0\\.", ".", x)
}

plot_predictions <- function(ideal_model_predictions, benchmarks_eval,
                             path, width, height) {
  get_ideal_pred <- function(x) {
    stopifnot(nrow(x) == 1)
    x$predictions[[1]] %>% 
      transmute(model = information_content,
                listener = surprisal_rating_mean)
  }
  
  dat <- list()
  
  dat$full <- ideal_model_predictions %>% 
    filter(ltm & ppm_order == 5 & viewpoint_group == "all") %>%
    get_ideal_pred() %>% 
    mutate(label = "full")
  
  dat$simple <- ideal_model_predictions %>% 
    filter(!ltm & ppm_order == 0 & viewpoint_group == "pc_chord_only") %>% 
    get_ideal_pred() %>% 
    mutate(label = "simple")
  
  dat$combined_benchmark <- benchmarks_eval %>% attr("combined_benchmark_data") %>% 
    mutate(label = "benchmark")
  
  dat_all <- bind_rows(dat) %>% 
    mutate(label = recode_factor(label,
                                 full = "Full ideal observer",
                                 benchmark = "Benchmark regression",
                                 simple = "Simplified ideal observer"))
  
  dat_corr <-
    dat_all %>% 
    group_by(label) %>% 
    summarise(
      spearman = cor(model, listener, method = "spearman") %>% 
        sprintf("%.2f", .) %>% drop_l0()
    )

  p <- dat_all %>% 
    ggplot(aes(model, listener)) +
    geom_bin2d() +
    scale_fill_viridis_c("Count") +
    scale_x_continuous("Model output") + 
    scale_y_continuous("Listener surprisal") +
    geom_text(data = dat_corr,
              aes(label = paste0("~italic(r[s]) == '", spearman, "'")),
              parse = TRUE,
              x = Inf, y = -Inf,
              hjust = 1.25, vjust = - 0.5) +
    facet_wrap(~ label, nrow = 1, scales = "free") +
    theme(aspect.ratio = 1)

  ggsave(filename = path, plot = p, width = width, height = height)
  p
}
         
cor_att <- function(x, y, x_se, y_se) {
  checkmate::qassert(x_se, "N1")
  checkmate::qassert(y_se, "N1")
  cor(x, y) / sqrt(sep_index(x, x_se) * sep_index(y, y_se))
}

sep_index <- function(x, x_se) {
  (var(x) - x_se ^ 2) / var(x)
}

