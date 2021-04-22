# d-social-kld

This repository contains replication code for D-SOCIAL-KLD scores ([Carroll, Primo, & Richter 2016](https://onlinelibrary.wiley.com/doi/full/10.1002/smj.2463)) using KLD data that runs through the year 2018.

## Files

Files in this directory include:

  - `README.md`: this file
  - `LICENSE`: license information
  - `01-cleanKLD.R`: a script for converting off-the-shelf KLD data into workable form for the `MCMCpack` routine; 
  - `02-getIdeals.R`: the script calling `MCMCpack::MCMCdynamicIRT1d`, which actually gets more than just the ideal points; 
  - `03-pullThetas.R`: the script pulling/cleaning the D-SOCIAL-KLD scores; 
  - `fn-get_start_values.R`: a helper function for starting values in the `MCMCpack::MCMCdynamicIRT1d` call; and
  - `d-social-kld_summary-stats.csv`: summary statistics for the posterior distributions.

The data with all 2500 pulls from the posterior are [too large to be added here](https://docs.github.com/en/github/managing-large-files/conditions-for-large-files), so please check [socialscores.org](http://socialscores.org/) for updates.

## Installation

The most recent analysis was conducted in April 2021.
The analysis was conducted in [R (Version 4.0.5)](https://cran.r-project.org/).
The only package dependencies are:

  - `tidyverse` (Version 1.3.0);
  - `foreach` (Version 1.5.1); and
  - `MCMCpack` (Version 1.5-0).
