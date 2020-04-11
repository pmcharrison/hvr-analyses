# Source this script to load functions and perform analyses.

# Preparation ####

library(tidyverse)
source("src/pop-behav/functions/ideal-observer-plots.R")

parent_dir <- "output/pop-behav"
out_dir <- file.path(parent_dir, "4-ideal-observer-plots")
R.utils::mkdirs(out_dir)
about <- readRDS(file.path(parent_dir, "about.rds"))
runs <- about$runs$ideal_regression %>% 
  mutate(results = map(ideal_regression_dir,
                       ~  readRDS(file.path(parent_dir, ., "results.rds"))))


# Parameters ####

chosen_order <- 5

viewpoint_labels <- runs$results[[1]]$viewpoint_labels %>% 
  mutate(viewpoint_label = if_else(viewpoint == "spec_sim_3", 
                                   "Spectral similarity", 
                                   viewpoint_label))

theme_set(theme_classic() + 
            theme(axis.text = element_text(colour = "black"),
                  axis.ticks = element_line(colour = "black")))

# Body ####

analyse_ppm_order(runs, out_dir)
analyse_feature_groups(runs, chosen_order, out_dir)
analyse_models(runs, chosen_order, viewpoint_labels, out_dir)
