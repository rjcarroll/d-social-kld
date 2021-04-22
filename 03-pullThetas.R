#
#
#    name: 03-pullThetas.R
#      by: rjc
#    date: 2021-04-23
#    what: pulls out the D-SOCIAL-KLD scores and makes them pretty.
#
#

#   ____________________________________________________________________________
#   preliminaries                                                           ####

    # clean up
    rm(list = ls())

    # packages
    library(tidyverse)
    library(MCMCpack)

    # read in the data
    load("data-MCMCout.RData")
    kld <- read_csv("firm-year-level.csv") %>% 
      dplyr::select(name, ticker, year)
    
    # clean up the workplace #
    rm(use, srs, theta.start, yrs)
    
#   ____________________________________________________________________________
#   pulling                                                                 ####
    
    th <- 
      as_tibble(out) %>% 
      dplyr::select(starts_with("theta")) %>% 
      rename_with(.fn = function(x) substring(x, first = 7))
    
    # summary stats
    # TODO: make this not so awful. good lord.
    th_summary <- 
      th %>% 
      pivot_longer(cols      = everything(),
                   names_to  = "key",
                   values_to = "value") %>%
      group_by(key) %>% 
      summarise(mean         = mean(value),
                sd           = sd(value),
                min          = min(value),
                percentile05 = quantile(value, probs = 0.05),
                percentile10 = quantile(value, probs = 0.10),
                percentile25 = quantile(value, probs = 0.25),
                percentile50 = quantile(value, probs = 0.50),
                percentile75 = quantile(value, probs = 0.75),
                percentile90 = quantile(value, probs = 0.90),
                percentile95 = quantile(value, probs = 0.95),
                max          = max(value))
    
    # clean up the key 
    th_summary <- 
      th_summary %>% 
      separate(col   = key,
               into  = c("ticker", "year"),
               sep   = "\\.t",
               extra = "merge") %>% 
      mutate(year = 1990 + as.integer(year)) %>% 
      arrange(ticker, year) %>% 
      left_join(y = kld,
                by = c("ticker", "year")) %>% 
      dplyr::select(name, ticker, year, everything())
    
    # write out the summary
    write_csv(x   = th_summary,
              file = "d-social-kld_summary-stats.csv")
    
#   end ## not run
      