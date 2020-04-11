library(tidyverse)
theme_set(theme_classic() +
            theme(axis.text = element_text(colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  strip.background = ggplot2::element_rect(
                    color = "#DEDEDE", size = 1,
                    fill = "#DEDEDE"
                  )))

parent_dir <- "output/pop-behav"

out_dir <- file.path(parent_dir, "6-continuous-descriptives")
R.utils::mkdirs(out_dir)

viewpoint_dir <- file.path(parent_dir, "0-viewpoints")
about <- readRDS(file.path(viewpoint_dir, "about.rds"))
mod <- readRDS(file.path(parent_dir, "3-ideal-regression", "run-1", "results.rds"))

viewpoint_labels <-
  mod$viewpoint_labels %>% 
  mutate(viewpoint_label = gsub(" \\(HL = [0-9]*\\)", "", viewpoint_label))

get_viewpoint_data <- function(viewpoint_dir, viewpoints, target_chord = 6L) {
  about <- readRDS(file.path(viewpoint_dir, "about.rds"))
  symbols <- readRDS(file.path(viewpoint_dir, "corpus.rds"))[about$seq_test] %>% 
    map_int(~ .[[target_chord]])
  
  stopifnot(all(about$seq_test == seq_along(about$seq_test)),
            length(viewpoints) > 1)
  
  all <- plyr::llply(about$seq_test, function(i) {
    x <- readRDS(file.path(viewpoint_dir, "viewpoints-test", paste0(i, ".rds")))
    x$continuous[viewpoints, target_chord, ]
  }, .progress = "text")
  
  map2(seq_along(all), symbols, function(seq_id, symbol) {
    event_mat <- all[[seq_id]]
    t(event_mat) %>% 
      as_tibble() %>% 
      add_column(pos = target_chord, .before = 1) %>% 
      add_column(symbol = seq_len(ncol(event_mat)), .before = 1) %>% 
      add_column(seq_id = seq_id, .before = 1) %>% 
      mutate(observed = symbol == !!symbol)
  }) %>% bind_rows()
}

viewpoints <- about$continuous_viewpoints %>% Filter(function(x) x != "num_pcs", .)

df <- get_viewpoint_data(viewpoint_dir, viewpoints)

set.seed(1)
p <- df %>%
  mutate(pi_dist = pmax(0, pi_dist + rnorm(length(pi_dist)))) %>%
  mutate(group = if_else(observed, "Observed", "Not observed") %>% 
           factor(levels = c("Not observed", "Observed"))) %>% 
  # select(- num_pcs) %>% 
  gather(key = "viewpoint", value = "value", viewpoints) %>% 
  mutate(viewpoint = plyr::mapvalues(viewpoint, 
                                     from = viewpoint_labels$viewpoint,
                                     to = viewpoint_labels$viewpoint_label,
                                     warn_missing = FALSE)) %>% 
  # slice(1:60000) %>%
  ggplot(aes(x = value, fill = group)) + 
  scale_x_continuous("Feature value") +
  scale_y_continuous("Density") +
  scale_fill_viridis_d(NULL) +
  # scale_fill_manual(NULL, values = c("#6ed1ff", "#B50000")) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ viewpoint, scales = "free") +
  theme(aspect.ratio = 1,
        legend.position = "bottom")
ggsave(file.path(out_dir, "continuous-descriptives.pdf"), plot = p, width = 6, height = 6)

df[, viewpoints] %>% 
  set_names(., plyr::mapvalues(viewpoints, 
                               from = viewpoint_labels$viewpoint,
                               to = viewpoint_labels$viewpoint_label,
                               warn_missing = FALSE)) %>% 
  cor() %>% 
  as.data.frame() %>% 
  write.csv(file.path(out_dir, "continuous-cor.csv"))
