library(hrep)
library(magrittr)
library(tidyverse)
loadNamespace("egg")

source("src/pop-behav/functions/plot-discrete-viewpoint.R")

ggplot2::theme_set(ggplot2::theme_classic() +
                     ggplot2::theme(axis.text = ggplot2::element_text(colour = "black"),
                                    axis.ticks = ggplot2::element_line(colour = "black")))

out_dir <- "output/pop-behav"
about <- readRDS(file.path(out_dir, "about.rds"))

plot_dir <- file.path(out_dir, "8-categorical-descriptives")
R.utils::mkdirs(plot_dir)

viewpoint_sequences <- readRDS(
  file.path(out_dir, "0-viewpoints", "viewpoints-training.rds")
)[about$indices$train]

plot_root_pc <- plot_discrete_viewpoint(viewpoint_sequences, "root_pc",
                                        alphabet_size = 12L,
                                        zero_indexed = TRUE)
plot_root_int <- plot_discrete_viewpoint(viewpoint_sequences, "root_int", 
                                         alphabet_size = 12L, 
                                         zero_indexed = TRUE)

plot_bass_pc <- plot_discrete_viewpoint(viewpoint_sequences, "bass_pc",
                                        alphabet_size = 12L,
                                        zero_indexed = TRUE)
plot_bass_int <- plot_discrete_viewpoint(viewpoint_sequences, "bass_int", 
                                         alphabet_size = 12L, 
                                         zero_indexed = TRUE)

combine_pair <- function(a, a_label, b, b_label, path) {
  p <- cowplot::plot_grid(
    a, b, 
    nrow = 1,
    labels = c(a_label, b_label),
    vjust = 1.5,
    scale = 0.85
  ) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(2, 0.5, 0.5, 0.5), "mm"))
  ggplot2::ggsave(path, plot = p, width = 8.55, height = 7)
  p
}

combine_pair(plot_root_pc, "A", 
             plot_root_int, "B",
             file.path(plot_dir, "roots.pdf"))

combine_pair(plot_bass_pc, "A", 
             plot_bass_int, "B",
             file.path(plot_dir, "bass.pdf"))


# Correlations
info <- about$runs$model_matrix %>% 
  filter(ltm & ppm_order == 5 & viewpoint_group == "categorical_only")

continuation_matrices <- readRDS(file.path(out_dir,
                                           info$model_matrix_dir, 
                                           "continuation-matrices.rds"))
predictors <- readRDS(file.path(out_dir,
                                info$model_matrix_dir,
                                "predictors.rds"))

vars <- predictors %>% 
  filter(discrete & class == "ltm") %>% 
  transmute(col = label,
            viewpoint = viewpoint,
            english = recode(viewpoint, 
                             bass_int = "Bass interval",
                             pc_set_rel_bass = "PC set rel. bass",
                             pc_set_rel_root = "PC set rel. root",
                             pc_chord = "PC chord",
                             pc_chord_rel_prev_bass = "PC chord rel. prev. bass",
                             pc_set = "PC set",
                             pc_set_rel_prev_bass = "PC set rel. prev. bass",
                             bass_pc_rel_root = "Bass PC rel. root",
                             root_int = "Root interval",
                             root_pc = "Root PC",
                             bass_pc = "Bass PC"
            ))
df <- map_dfr(continuation_matrices, ~ as_tibble(.[, vars$col]))
cor_mat <- cor(df)
dist_mat <- as.dist(1 - cor_mat)

pheatmap::pheatmap(cor_mat,
                   color = viridis::viridis(100),
                   border_color = "black",
                   legend_breaks = seq(from = 0, to = 1, by = 0.2),
                   cellwidth = 20, cellheight = 20,
                   clustering_distance_rows = dist_mat, # "euclidean",
                   clustering_distance_cols = dist_mat, # "euclidean",
                   labels_row = vars$english,
                   labels_col = vars$english,
                   filename = file.path(plot_dir, "cor-mat.pdf"),
                   width = 6, height = 6, dpi = 300
)

