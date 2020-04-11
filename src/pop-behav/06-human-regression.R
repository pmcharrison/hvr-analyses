# library(tidyverse)
# library(hrep)
# library(hvr)
# source("src/pop-behav/functions/human-regression.R")
# 
# parent_dir <- "output/pop-behav"
# out_dir <- file.path(parent_dir, "5-human-regression")
# chosen_order <- 4
# 
# run_human_regressions(parent_dir,
#                       out_dir,
#                       chosen_order = chosen_order,
#                       perm_int_reps = 1,
#                       max_eval = 1000L) # 1000L
