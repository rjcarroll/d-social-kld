# d-social-kld

This repository contains replication code for D-SOCIAL-KLD scores ([Carroll, Primo, & Richter 2016](https://onlinelibrary.wiley.com/doi/full/10.1002/smj.2463)) using KLD data that runs through the year 2018.

## Files

Files in this directory include:

  - `README.md`: this file;
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

## Notes

  - There are well-known and oft-lamented classes between `dplyr::select` and the `MASS` package called in `MCMCpack`. We've tried to minimize errors, but the occasional sequencing error may have slipped through here or there.
  - Several firms from distinct KLD "domiciles" share the same ticker symbol. In these cases, the newer-to-the-dataset firm is assigned a random four-character tag placed after a double-tilde separator. So, if two firms both have ticker symbol `Z`, then the older-to-the-dataset firm remains `Z` and the newer-to-the-dataset firm is assigned ticker `Z~~abcd`, where `abcd` may be any random tag of four lowercase letters. Those that read the code in `01-cleanKLD.R` carefully will see that several precautions were taken to ensure continuity across these tags. However, should you find any errors, please do not hesitate to drop a line or issue a pull request. 
