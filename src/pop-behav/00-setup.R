library(tidyverse)
library(testthat)
library(hrep)

source("src/pop-behav/functions/about.R")

out_dir <- "output/pop-behav"
ppm_orders <- 0:5
ppm_training_sizes <- 439L # maximum = 439
poly_degrees <- 4L

viewpoint_group <- c("all", "continuous_only", "categorical_only")

df <- get_stimulus_data('input/pop-behav/pop-behav.csv')
stimulus_ids <- sort(unique(df$stimulus_id))
test_corpus <- get_test_corpus(stimulus_ids, df)
train_corpus <- get_train_corpus(df)
combined_corpus <- c(test_corpus, train_corpus)
test_ids <- seq_along(test_corpus)
train_ids <- seq(from = length(test_corpus) + 1L,
                 to = length(combined_corpus))
runs <- get_runs(ppm_orders,
                 ppm_training_sizes,
                 poly_degrees,
                 train_ids,
                 viewpoint_group)

save_about(out_dir, df, test_corpus, train_corpus, combined_corpus, 
           test_ids, train_ids, stimulus_ids, runs)
