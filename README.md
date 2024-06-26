# hvr-analyses

<!-- badges: start -->
<!-- badges: end -->

This repository contains source code for the harmonic expectation chapter
of Peter M. C. Harrison's (2020) PhD thesis. 

The following code installs the custom R packages required for these
analyses:

``` r
install.packages("remotes")
remotes::install_github(c("pmcharrison/hrep",
                          "pmcharrison/incon",
                          "pmcharrison/hvr"))
```

Other packages can be installed as needed using `install.packages`.

The code is intended to be run by sourcing the files in `src/pop-behav`
in ascending numeric sequence. They save their outputs to `output`.
Each of these files can be run from a blank workspace,
assuming that the previous files have been run already
and their results saved to `output`.

The original code for `04-ideal-regression.R` computed permutation-based 
feature importance metrics that were very time-consuming to evaluate.
These have been expedited by default by setting `perm_int_reps = 1` in
`04-ideal-regression.R`. To compute these properly one would have to increase
that number to e.g. 100.
