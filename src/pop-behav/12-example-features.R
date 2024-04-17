library(tidyverse)
library(hrep)
theme_set(theme_classic() +
            theme(axis.text = element_text(colour = "black"),
                  axis.ticks = element_line(colour = "black"),
                  strip.background = ggplot2::element_rect(
                    color = "#DEDEDE", size = 1,
                    fill = "#DEDEDE"
                  )))

parent_dir <- "output/pop-behav"

out_dir <- file.path(parent_dir, "9-example-features")
R.utils::mkdirs(out_dir)

viewpoint_dir <- file.path(parent_dir, "0-viewpoints")
about <- readRDS(file.path(viewpoint_dir, "about.rds"))
corpus <- readRDS(file.path(viewpoint_dir, "corpus.rds"))

viewpoint_labels <- about$viewpoint_labels

dat <- readRDS(file.path(viewpoint_dir, "viewpoints-training.rds"))

x <- dat[[1]]

x

format_set <- function(x) {
  paste0("\\{", paste(as.numeric(x), collapse = ", "), "\\}")
}

format_pc_set <- function(x) {
  if (is.na(x)) "NA" else
    decode(x, "pc_set")[[1]] %>%
    format_set()
}

format_pc_chord <- function(x) {
  if (is.na(x)) "NA" else
    decode(x, "pc_chord")[[1]] %>%
    {
      sprintf("(%i, \\{%s\\})", 
              get_bass_pc(.), 
              paste(get_non_bass_pc(.), collapse = ", "))
    }
}

format_pc <- function(x) {
  if (is.na(x)) "NA" else as.character(x - 1)
}

format_continuous <- function(x) {
  if (is.na(x)) "NA" else sprintf("%.2f", x)
}

format_integer <- function(x) {
  if (is.na(x)) "NA" else sprintf("%i", round(x))
}

format_pc_set_rel_bass <- function(x) {
  decode(x, "pc_chord_type")[[1]] %>% 
    format_set()
}

format_pc_set_rel_root <- function(x) {
  mapping <- tibble(
    pc_set_id = as.integer(hvrmap::map_pc_chord$pc_set_rel_root_id %>% levels()),
    pc_set_rel_root_id = seq_along(pc_set_id)
  )
  pc_set_id <- mapping$pc_set_id[x]
  decode(pc_set_id, "pc_set")[[1]] %>% 
    format_set()
}


formatters <- list(
  bass_int = format_pc,
  bass_pc = format_pc,
  bass_pc_rel_root = format_pc,
  pc_chord = format_pc_chord,
  pc_chord_rel_prev_bass = format_pc_chord,
  pc_set = format_pc_set,
  pc_set_rel_bass = format_pc_set_rel_bass,
  pc_set_rel_prev_bass = format_pc_set,
  pc_set_rel_root = format_pc_set_rel_root,
  root_int = format_pc,
  root_pc = format_pc,
  har_18_harmonicity = format_continuous,
  hutch_78_roughness = format_continuous,
  num_pcs = format_integer,
  pi_dist = format_integer,
  spec_sim_3 = format_continuous
)

data <- x$discrete

format_features <- function(data) {
  features <- rownames(data)
  feature_labels <- plyr::mapvalues(features, 
                                    from = viewpoint_labels$viewpoint,
                                    to = viewpoint_labels$viewpoint_label,
                                    warn_missing = FALSE)
  
  map(seq_len(nrow(data)), function(i) {
    map_chr(data[i, ], 
            formatters[[features[i]]]) %>% 
      as.list() %>% 
      {set_names(., paste("", seq_along(.)))} %>%  
      as_tibble()
  }) %>% 
    bind_rows() %>% 
    add_column(Feature = feature_labels, .before = 1)
}

stim_id <- 5


list(
  description = corpus[[stim_id]] %>% metadata() %>% `$`(description),
  
  # Originally this code was written to visualise continuous features too,
  # but for some reason these are rendering as NULL, so now we just 
  # subset by the discrete features below.
  features = map(dat[[stim_id]]["discrete"], format_features)
) %>% 
  saveRDS(file.path(out_dir, "example-features.rds"))

