
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RankspeQ

<!-- badges: start -->
<!-- badges: end -->

RankspeQ is a novel R package developed to evaluate genotype performance
and support selection-driven decisions based on leaf traits and also
environment-related variables measured by the MultispeQ device. The
presented software consists of 3 main functions: the data cleaning,
genotype ranking based on its performance and comparison of accessions
against grain yield or another crop trait. Optionally, the evaluation
can be done by groups which can be defined by the user such as
genepools, families, races, etc. The software development as well as the
data evaluation was made with datasets of Phaseolus spp. experiments.
However, R code - with easy modifications - can be used on any other
type of crop. The tool will help to understand the hidden potential of
MultispeQ equipment and identify important crop traits useful in a
genotype characterization in particular environments. The tool has
direct potential for physiologists, breeders etc. as he identifies the
best performing accessions, however, also target false positive results
with high yield but low photosynthetic performance. We also propose to
use a new efficiency index (Phi Index). Further updates will include new
algorithms (e.g. trait heritability) and a shiny interface to make the
software user friendly.

## Installation

You can install the development version of RankspeQ from
[GitHub](https://github.com/jssotob/RankspeQ) with:

``` r
# install.packages("devtools")
devtools::install_github("jssotob/RankspeQ")
```
