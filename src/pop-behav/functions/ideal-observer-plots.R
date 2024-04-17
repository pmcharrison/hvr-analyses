library(tidyverse)
library(magrittr)

analyse_ppm_order <- function(runs, out_dir) {
  runs %>% 
    filter(viewpoint_group == "all") %>% 
    transmute(ppm_order = ppm_order, 
              ltm = ltm,
              cross_entropy_bits = map_dbl(results, "cost")) %T>% 
    write_csv(file.path(out_dir, "ppm-order.csv"))
}

analyse_feature_groups <- function(runs, chosen_order, out_dir) {
  runs %>% 
    filter(ppm_order == chosen_order |
             viewpoint_group == "continuous_only") %>%
    transmute(ppm_order = ppm_order,
              feature_group = viewpoint_group,
              ltm = ltm,
              cross_entropy_bits = map_dbl(results, "cost"),
              aic = map_dbl(results, "aic")) %T>% 
    write_csv(file.path(out_dir, "feature-groups.csv"))
}

analyse_models <- function(runs, chosen_order, viewpoint_labels, out_dir) {
  model_plots <- get_model_plots(runs, chosen_order, viewpoint_labels)
  export_marginals(model_plots, out_dir)
  export_perm_int(model_plots, out_dir)
  export_categorical_weights(runs, out_dir)
}

get_model_plots <- function(runs, chosen_order, viewpoint_labels) {
  c("all", "continuous_only") %>% 
    setNames(., .) %>% 
    map(function(x) {
      mod <- runs %>% 
        filter(viewpoint_group == x) %>% 
        filter(if_else(viewpoint_group == "continuous_only", 
                       TRUE,
                       ltm & ppm_order == chosen_order)) %>% 
        pull(results) %>% 
        `[[`(1)
      list(
        marginals = hvr::plot_marginals(mod, viewpoint_labels = viewpoint_labels, fill = "#21908CFF") +
          ggplot2::theme(strip.background = ggplot2::element_rect(
            color = "#DEDEDE", linewidth = 1,
            fill = "#DEDEDE"
          )),
        perm_int = hvr::plot_perm_int(mod, 
                                      labels = viewpoint_labels,
                                      discrete_label = "High-level",
                                      continuous_label = "Low-level") + 
          theme(legend.position = if (x == "all") c(0.8, 0.85) else "none")
      ) %>% 
        map(`+`, theme(aspect.ratio = 1))
    })
}

export_marginals <- function(model_plots, out_dir) {
  export_marginal_plot(model_plots, out_dir)
  export_marginal_data(model_plots, out_dir)
}

export_marginal_plot <- function(model_plots, out_dir) {
  map(model_plots, "marginals") %>% 
    cowplot::plot_grid(plotlist = ., nrow = 1, labels = "AUTO", scale = 0.9) %>% 
    ggsave(plot = ., file.path(out_dir, "marginals.pdf"), width = 8.5, height = 4.75)
}

export_marginal_data <- function(model_plots, out_dir) {
  map(model_plots, "marginals") %>% 
    map(~ attr(., "data")) %>% 
    map2(names(model_plots), 
         ~ add_column(.x, model_type = .y, .before = 1)) %>% 
    bind_rows() %>% 
    rename(feature = viewpoint, 
           feature_label = viewpoint_label) %>% 
    write_csv(file.path(out_dir, "marginals.csv"))
}

export_perm_int <- function(model_plots, out_dir) {
  export_perm_int_plot(model_plots, out_dir)
  export_perm_int_data(model_plots, out_dir)
}

export_perm_int_plot <- function(model_plots, out_dir) {
  map(model_plots, "perm_int") %>% 
    cowplot::plot_grid(plotlist = ., nrow = 1, labels = "AUTO", scale = 1) %>% 
    ggsave(plot = ., file.path(out_dir, "perm-int.png"), width = 9, height = 4.1, dpi = 300)
}

export_perm_int_data <- function(model_plots, out_dir) { 
  map(model_plots, "perm_int") %>% 
    map(~ attr(., "data")) %>% 
    map2(names(model_plots), 
         ~ add_column(.x, model_type = .y, .before = 1)) %>% 
    map(~ mutate(., viewpoint = as.character(viewpoint))) %>% 
    bind_rows() %>% 
    rename(feature = viewpoint) %>% 
    write_csv(file.path(out_dir, "perm-int.csv"))
}

export_categorical_weights <- function(runs, out_dir) {
  mod <- runs %>% 
    filter(ppm_order == chosen_order & viewpoint_group == "all") %>% 
    pull(results) %>% 
    `[[`(1)
  
  p <- hvr::plot_discrete_weights(mod,
                                  x_lab = "Feature", 
                                  colours = viridis::viridis(2, direction = -1)) + 
    theme(aspect.ratio = 1)
  data <- attr(p, "data") %>% 
    rename(categorical = discrete,
           feature = viewpoint,
           feature_label = viewpoint_label)
  
  ggsave(plot = p, file.path(out_dir, "categorical-weights.pdf"), width = 6, height = 3.5)
  write_csv(data, file.path(out_dir, "categorical-weights.csv"))
}
