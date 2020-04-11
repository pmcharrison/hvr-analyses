library(tidyverse)
library(hvr)
library(hrep)

out_dir <- "output/pop-behav/10-example-seq-analyses"

model_dataset(
  corpus_test = hcorp::popular_1[1:5],
  corpus_pretrain = hcorp::popular_1[seq(to = 739, length.out = 700)],
  output_dir = out_dir,
  weights = vr_config$hvr_1$weights,
  viewpoints = vr_config$hvr_1$viewpoints,
  poly_degree = 4
)

if (FALSE) {
  # Inspecting results
  x <- readRDS(file.path(out_dir, "4-predictions", "predictions.rds")) %>% 
    group_by(seq_id) %>% 
    mutate(voiced = voicer::voice(chord, 
                                  opt = voicer::voice_opt(min_octave = -1,
                                                          max_octave = 0,
                                                          max_notes = 6))) %>% 
    group_split()
  
  x[[1]] %>% 
    {hrep::view(.$voiced, 
                chords_per_line = 4,
                annotate = sprintf("%.1f", .$information_content))}
}


