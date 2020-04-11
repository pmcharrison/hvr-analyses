# all.rds

This RDS file provides detailed results for the behavioural experiment. It also includes some preliminary modelling analyses, but these are incomplete and should be ignored.

# pop-behav.csv

## Summary

This dataset provides chord-level features for our 300 stimuli, including aggregated surprisal judgements from participants. The dataset was sourced from the HarmonyExpectancyEEG project.

## Columns

`stimulus-id` (integer) - Unique identifier for each chord sequence. Ranges from 521 to 580.

`stimulus_label` (string) - Longer unique identifier for each chord sequence (e.g. *id=521_genre=popular_c-id=508_e-id=30_ic_category_1*).

`title` (string) - Artist and song title (e.g. *Creedence Clearwater Revival: Bad Moon Rising*).

`dataset_id` (integer) - ID of the musical corpus in Peter's IDyOM implementation. Should be equal to 2 for all pieces in this popular music corpus.

`composition_id` (integer) - ID of the song from which the chord sequnce is excerpted in Peter's IDyOM implementation (0-indexed).

`chord_id` (integer) - Position of the current chord within the chord sequence (1-indexed), ranges from 1 to 8.

`midi` (integer vector) - MIDI note numbers present in the current chord as played to the participant, e.g. '55 62 71'.

`pc_set` (integer vector) - Pitch-class set of the current chord, e.g. `0 4 7`, where `0` corresponds to 'C'.

`surprisal_rating_mean` (numeric) - Mean normalised surprisal rating as given by participants in the preceding behavioural study (Harrison & Pearce, ICMPC 2018). This data is only available for the 6th chord of each chord sequence.

`surprisal_rating_sd` (numeric) - Standard deviation of normalised surprisal ratings across participants (see `surprisal_rating_mean`). 

`surprisal_rating_n` (integer) - Number of participant ratings contributing to this chord's value of `surprisal_rating_mean`.

`surprisal_rating_se` (numeric) - Standard error of `surprisal_rating_mean` as estimated using the central limit theorem.

`probability` (numeric) - Probability of the chord as estimated by Harrison & Pearce's viewpoint regression model.

`information_content` (numeric) - Information content of the chord as estimated by Harrison & Pearce's viewpoint regression model; computed as the negative base-2 logarithm of `probability`.

`milne_2011` (numeric) - Spectral distance between the current chord and the preceding chord as computed by Milne et al. (2011)'s pitch-class model (see Harrison & Pearce, ICMPC 2018).

`onset` (numeric) - Temporal onset of the current chord (seconds), relative to the beginning of the audio file.

`onset` (numeric) - Temporal offset of the current chord (seconds), relative to the beginning of the audio file.

`leman_2000_v1` (numeric) - Tonal contextuality of the current chord as computed by Leman's (2000) Periodicity Pitch model and converted to a z-score by Fisher's z-transformation after Bigand et al. (2014). This dataset contains Leman model outputs for a variety of tuning parameters for the local and global decay constants, after Bigand et al. (2014). For up-to-date definitions of these tuning parameters, see `output/03-astm/astm-features.rds`, which can be opened in R. This column corresponds to a local decay parameter of 0.1s and a global decay parameter of 1.5 s.

`leman_2000_v2` (numeric) - `leman_2000_v1` but with local decay of 0.1 s and global decay of 2.5 s.

`leman_2000_v3` (numeric) - `leman_2000_v1` but with local decay of 0.1 s and global decay of 4.0 s.


`leman_2000_v4` (numeric) - `leman_2000_v1` but with local decay of 0.5 s and global decay of 1.5 s.


`leman_2000_v5` (numeric) - `leman_2000_v1` but with local decay of 0.5 s and global decay of 2.5 s.

`leman_2000_v6` (numeric) - `leman_2000_v1` but with local decay of 0.5 s and global decay of 4.0 s.

`collins_2014` (numeric) - output of the model from Collins et al. (2014) as implemented in the Janata Lab Music Toolbox and Peter's [matlab-expectation-models](https://github.com/pmcharrison/matlab-expectation-models) repository.

`melody_midi` (numeric) - highest MIDI note in the chord, termed the melody note

`melody_information_content` (numeric) - information content (i.e. surprisal) of the chord's melody note as determined by IDyOM (higher values mean higher surprise) 

`melody_entropy` (numeric) - entropy (i.e. predictive uncertainty) of IDyOM when trying to predict the melody note of this chord

`melody_dist` (numeric) - (unsigned) distance between this chord's melody note and the preceding chord's melody note, in semitones.

`vl_dist` (numeric) - minimal voice-leading distance between this chord and the preceding chord as estimated by Tymoczko's (2006) algorithm. Takes NA values for the first chord in each sequence. High distances mean that pitches have to move far to connect successive chords.

`incon` (numeric) - consonance estimate for the chord as produced by Harrison & Pearce's (2019) composite consonance algorithm (higher values indicate greater consonance).

## References

Bigand, E., Delbé, C., Poulin-Charronnat, B., Leman, M., & Tillmann, B. (2014). Empirical evidence for musical syntax processing? Computer simulations reveal the contribution of auditory short-term memory. Frontiers in Systems Neuroscience, 8. https://doi.org/10.3389/fnsys.2014.00094

Collins, T., Tillmann, B., Barrett, F. S., Delbé, C., & Janata, P. (2014). A combined model of sensory and cognitive representations underlying tonal expectations in music: From audio signals to behavior. Psychological Review, 121(1), 33–65. https://doi.org/10.1037/a0034695

Harrison, P. M. C., & Pearce, M. T. (2018). Dissociating sensory and cognitive theories of harmony perception through computational modeling. In Proceedings of ICMPC15/ESCOM10. Graz, Austria. https://doi.org/10.31234/osf.io/wgjyv

Harrison, P. M. C., & Pearce, M. T. (2019). Instantaneous consonance in the perception and composition of Western music. PsyArXiv. https://doi.org/10.31234/osf.io/6jsug

Leman, M. (2000). An auditory model of the role of short-term memory in probe-tone ratings. Music Perception, 17(4), 481–509.

Milne, A. J., Sethares, W. A., Laney, R., & Sharp, D. B. (2011). Modelling the similarity of pitch collections with expectation tensors. Journal of Mathematics and Music, 5(1), 1–20. https://doi.org/10.1080/17459737.2011.573678

Tymoczko, D. (2006). The geometry of musical chords. Science, 313(5783), 72–74. https://doi.org/10.1126/science.1126287
