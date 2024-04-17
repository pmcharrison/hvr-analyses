library(tidyverse)
library(memoise)
library(magrittr)
source("src/pop-behav/functions/human-regression.R")
source("src/pop-behav/functions/eval-human-models.R")

parent_dir <- "output/pop-behav"
out_dir <- file.path(parent_dir, "5-eval-human-models")
R.utils::mkdirs(out_dir)
about <- readRDS(file.path(parent_dir, "about.rds"))
cache_dir <- file.path(parent_dir, "cache")
R.utils::mkdirs(cache_dir)

get_ideal_model_predictions <- memoise(get_ideal_model_predictions, 
                                       cache = cache_filesystem(path = cache_dir))
get_benchmarks <- memoise(get_benchmarks,
                          cache = cache_filesystem(path = cache_dir))

theme_set(theme_classic() +
            theme(axis.text = element_text(colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  strip.background = ggplot2::element_rect(
                    color = "#DEDEDE", size = 1,
                    fill = "#DEDEDE"
                  )))

###

ideal_model_predictions <- get_ideal_model_predictions(parent_dir,
                                                       about,
                                                       max_information_content = 15) %>% 
  save_ideal_model_predictions(out_dir)

benchmarks <- get_benchmarks() %T>% 
  write_csv(file.path(out_dir, "benchmarks.csv"))

benchmarks_eval <- eval_benchmarks(benchmarks) %T>% 
  write_csv(file.path(out_dir, "benchmarks-eval.csv"))

eval_ideal_vs_benchmarks <- get_eval_ideal_vs_benchmarks(benchmarks_eval, 
                                                         ideal_model_predictions) %T>%
  write_csv(file.path(out_dir, "eval-ideal-vs-benchmarks.csv"))

plot_predictions(ideal_model_predictions, benchmarks_eval, 
                 file.path(out_dir, "predictions.pdf"),
                 width = 7.5, height = 3)

ids <- list(
  full = ideal_model_predictions %>%
    filter(ltm & ppm_order == 5 & viewpoint_group == "all") %>%
    pull(ideal_regression_run_id),
  
  simple = ideal_model_predictions %>%
    filter(!ltm & ppm_order == 0 & viewpoint_group == "pc_chord_only") %>%
    pull(ideal_regression_run_id)
)

cor_compare <- get_cor_compare(ideal_model_predictions) %T>%
  write_csv(file.path(out_dir, "cor-compare.csv"))

cor_compare %>% 
  filter(i == ids$full & j == ids$simple)

cor_models(ids$full, ids$simple, ideal_model_predictions, method = "pearson")

regress_models(ids$full, ids$simple, ideal_model_predictions) %>%
  {list(
    mod = broom::glance(.),
    pred = broom::tidy(., conf.int = TRUE)
  )}

# 
# if (FALSE) {
#   ideal_model_predictions %>% 
#     select(ideal_regression_run_id, viewpoint_group, ltm,
#            ppm_order, spearman_rho, pearson_r, cross_entropy) %>% 
#     View()
# }
# 
# cor_model(6, ideal_model_predictions)
# cor_model(1, ideal_model_predictions)
# cor_model(19, ideal_model_predictions)
# cor_model(7, ideal_model_predictions)
# cor_model(26, ideal_model_predictions)
# cor_model(32, ideal_model_predictions)
# 
# benchmark_model(6, ideal_model_predictions, "input/pop-behav/pop-behav.csv")
# 
# cor_models(6, 7, ideal_model_predictions)
# cor_models(19, 7, ideal_model_predictions)
# cor_models(6, 32, ideal_model_predictions)
# 
# plot_model(6, ideal_model_predictions)
# 
# regress_models(6, 32, ideal_model_predictions) %>% summary()
# regress_models(6, 19, ideal_model_predictions) %>% summary()
