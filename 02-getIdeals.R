#
#
#    name: 02-getIdeals.R
#      by: rjc
#    date: 2021-04-23
#    what: this runs the MCMC routine
#
#

#   ____________________________________________________________________________
#   introductories                                                          ####

    # clean up 
    rm(list = ls())

    # packages
    library(tidyverse)
    library(MCMCpack)

    # helper function for priors
    source("fn-get_start_values.R")

#   ____________________________________________________________________________
#   load in and prep data                                                   ####

    load("data-reshaped.RData")

    srs <- get_start_values(dat)

    dat            <- as.matrix(dat)
    row.names(dat) <- allTicks
    
    # priors
    theta.start <- srs$firmper * 15
    
    # anchors
    theta.start[which(rownames(dat) == "HAL")]  <- -3
    theta.start[which(rownames(dat) == "MSFT")] <- +3
    
#   ____________________________________________________________________________
#   run it                                                                  ####
    
    set.seed(61802)
    # make sure you've got a lot of memory clear and nothing all that pressing
    # to do in the next 24-96 hours.
    out <- MCMCdynamicIRT1d(datamatrix        = dat,
                            item.time.map     = yrs,
                            theta.start       = theta.start,
                            alpha.start       = srs$metper,
                            beta.start        = 0,
                            mcmc              = 5000,
                            burnin            = 1000,
                            thin              = 2,
                            verbose           = 1,
                            seed              = NA,
                            theta.constraints = list(HAL = "-", MSFT = "+"))
    # this will throw a warning on the is.na(alpha.start) criterion---i think
    # this isn't on the user end.

    save.image(file = "data-MCMCout.RData")

#   end ## not run