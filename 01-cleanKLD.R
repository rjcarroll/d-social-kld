#
#
#    name: 01-cleanKLD.R
#      by: rjc
#    date: 2021-04-22
#    what: this cleans and reshapes the KLD data
#
#

#   ____________________________________________________________________________
#   introductories                                                          ####

    # clean up
    rm(list = ls())

    # packages
    library(tidyverse)
    library(foreach)

#   ____________________________________________________________________________
#   load in and clean data                                                  ####

    # read in data
    kld <- 
      read_csv(file      = "KLD_stats_all_through_2018.csv",
               na        = c("", "R", "NA"),
               # the following variables are wonky and not of use:
               col_types = cols('CUSIP'            = col_skip(),
                                'legacy_companyID' = col_skip(),
                                'domicile'         = col_skip(),
                                'issuerid'         = col_skip())) %>% 
      # remaining parsing problems are all when HUM_str_num = 2, so:
      mutate(HUM_str_num = if_else(condition = HUM_str_num >= 1, 
                                   true      = 1,
                                   false     = 0)) %>%
      # get rid of cumulative totals:
      select(-contains("num")) %>% 
      # put things in order:
      select(CompanyName, Ticker, year, everything()) %>% 
      # rename to save some annoyances:
      rename('name'   = CompanyName,
             'ticker' = Ticker)
    
    # fix the logicals
    fixLogical <- function (x)
    {
      
      y <- rep(NA, length(x))
      y[x == FALSE] <- 0
      y[x == TRUE]  <- 1
      
      as.integer(y)
      
    }
    
    cleanTwo <- function (x)
    {
      
      y <- x
      y[x >= 2] <- 1
      
      as.integer(y)
      
    }
    
    kld <- 
      kld %>% 
      mutate(across(.cols = where(is_logical),
                    .fns  = fixLogical)) %>% 
      mutate(across(.cols = where(is_double),
                    .fns  = as.integer))
    
    kld <- 
      kld %>% 
      mutate(across(.cols = where(is_double) & !c(year),
                    .fns  = cleanTwo))
    
    # rather than parse through all the different domiciles, let's just 
    # get the tickers of firms from 1991-2012 data (ensures this is similar in
    # kind with original D-SOCIAL-KLD)
    universe <- 
      kld %>% 
      filter(year <= 2012) %>% 
      distinct(ticker) %>% 
      pull(ticker)
      
    kld <- 
      filter(kld, ticker %in% universe)
    
    rm(universe)
    
##  ............................................................................
##  some problematic tickers/names                                          ####
    
    # having USX as "X" throws off some statistics later on, so:
    kld <- 
      kld %>%
      mutate(ticker = ifelse(ticker == "X", "USX", ticker)) 
    
    # LLL and MKL both have bad duplicates
    kld$ticker[kld$name == "lululemon athletica inc."]     <- "LLLM"
    kld$ticker[kld$name == "Alterra Capital Holdings Ltd"] <- "MKLA"
    
    # mistake with Heritage Financial Corporation versus group
    kld$ticker[kld$name == "Heritage Financial Corporation"] <- "HFWA"
    
    # several missing tickers:
    # Ryder System is R, which created a problem from the read-in:
    kld$ticker[kld$name %in% c("Ryder System, Inc.",
                               "RYDER SYSTEM, INC.")] <- "R"
    # one missing ticker in AuthenTec, Inc.
    kld$ticker[kld$name == "AuthenTec, Inc."] <- "AUTH"
    # one missing ticker in Benihana Inc
    kld$ticker[kld$name == "Benihana Inc"] <- "BNHNA"
    # one missing ticker in MEDTOX Scientific, Inc.
    kld$ticker[kld$name == "MEDTOX Scientific, Inc."] <- "MTOX"
    # national bank of canada has ticker NA, which created a problem
    kld$ticker[kld$name %in% c("National Bank of Canada",
                               "NATIONAL BANK OF CANADA",
                               "BANQUE NATIONALE DU CANADA")] <- "NA"
    # no tickers for Romarco Minerals Inc.
    kld$ticker[kld$name == "Romarco Minerals Inc."] <- "RTRAF"
    # no tickers for AKER SOLUTIONS ASA 
    kld$ticker[kld$name == "AKER SOLUTIONS ASA"] <- "AKSO"
    # capital property fund is a problem: CPL is doubled, and some NAs.
    kld$ticker[kld$name %in% c("CAPITAL PROPERTY FUND",
                               "Capital Property Fund Ltd")] <- "CPL:SJ"
    # no tickers for Intra-Cellular Therapies Inc
    kld$ticker[kld$name == "Intra-Cellular Therapies Inc"] <- "ITCI"
    # no tickers for MGM GROWTH PROPERTIES LLC
    kld$ticker[kld$name == "MGM GROWTH PROPERTIES LLC"] <- "MGP"
    # problems with Allergan
    kld$ticker[kld$name %in% c("Allergan, Inc.", "Allergan plc")] <- "AGN"
    kld$ticker[kld$name == "AEGON N.V."] <- "AEG~~uixl"
    # problems with choice point
    kld$ticker[kld$name %in% c("Cyfrowy Polsat S.A.", 
                               "CYFROWY POLSAT SPOLKA AKCYJNA")] <- "CPS~~xuqn"
    kld$ticker[kld$name == "COOPER-STANDARD HOLDINGS INC."] <- "CPS~~iuvh"
    # problems with Ladenburg Thalmann
    kld$ticker[kld$name == "GRUPA LOTOS SPOLKA AKCYJNA"] <- "LTS~~ssil"
    
#   ____________________________________________________________________________
#   learn which metrics and firms are in which year                         ####
    
##  ............................................................................
##  which metrics are in which year?                                        ####
    
    has_data <- function (x) any(!is.na(x))
    all_distinct <- function (x, negate = FALSE) 
    {
      
      foo <- identical(sort(x), sort(unique(x)))
      
      if (negate == FALSE) foo 
      if (negate == TRUE) !foo
      
    }
      
    
    mets_by_year <- 
      foreach (i = min(kld$year):max(kld$year)) %do%
      {
        
        filter(kld, year == i) %>% 
          select(-name, -ticker, -year) %>% 
          select_if(has_data) %>% 
          names(.) 

      }
    rm(i)
    
    # check to see if all metrics are distinct within year
    lapply(mets_by_year, all_distinct) %>% 
      unlist() %>% 
      all()
    # they are indeed
    
##  ............................................................................
##  which tickers are in which year?                                        ####
    
    firms_in_year <- 
      foreach (i = min(kld$year):max(kld$year)) %do% 
      {
        
        filter(kld, year == i) %>% 
          pull(ticker) %>% 
          sort(.)
        
      }
    rm(i)
    
    # check to see if all tickers are distinct within year
    lapply(firms_in_year, all_distinct, negate = TRUE) %>% 
      unlist() %>% 
      which(.)
    # they are not :-(.
    
    # 10 is a bad year
    filter(kld, year == 1990 + 10) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num))
    # for some reason, we have four copies of the same row for RIG. we can just
    # unique our way out of that one.
    firms_in_year[[10]] <- sort(unique(firms_in_year[[10]]))
    
    # 15 is a bad year
    filter(kld, year == 1990 + 15) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num))
    # same with BWNG in year 15 we can unique our way out of it.
    firms_in_year[[15]] <- sort(unique(firms_in_year[[15]]))
    
    # 16 is a bad year
    filter(kld, year == 1990 + 16) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num))
    # same with H in year 16.
    firms_in_year[[16]] <- sort(unique(firms_in_year[[16]]))
    
    # 17 is a bad year
    filter(kld, year == 1990 + 17) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num))
    # same with IVZ_LN in year 17.
    firms_in_year[[17]] <- sort(unique(firms_in_year[[17]]))
    
    # 20 is a bad year
    filter(kld, year == 1990 + 20) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num))
    # same with BNHNA, CENT, HEI, and UBA in year 20 (checked all---duplicates)
    firms_in_year[[20]] <- sort(unique(firms_in_year[[20]]))
    
    # 22 is a bad year
    filter(kld, year == 1990 + 22) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num))
    # same with KCAP and UAM in year 22
    firms_in_year[[22]] <- sort(unique(firms_in_year[[22]]))
    
    # 23 is a bad year
    filter(kld, year == 1990 + 23) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num)) -> bad23
    # MANY bad tickers. 190.
    filter(kld, year == 1990 + 23 & ticker %in% bad23$ticker) %>% 
      group_by(ticker) %>% 
      summarise(numNames = length(unique(name)))
    # and each of these has multiple company names. it's fine to just index
    # these now (assuming we use a rule that will work for future years). that
    # said, we also need to ensure that previous tickers are linked to this one.
    for (i in 1:nrow(bad23))
    {
      
      # see if the ticker has been represented earlier.
      oldI <- filter(kld, year < 1990 + 23 & ticker == bad23$ticker[i]) %>%
        select(name, ticker, year)
      
      # figure out which of the duplicated tickers are new
      if (nrow(oldI) > 0)
      {
        
        newNamesI <- 
          filter(kld, year == 1990 + 23 & ticker == bad23$ticker[i]) %>% 
          distinct(name) %>% 
          pull(name)
        
        notInOld <- setdiff(newNamesI, oldI$name)
        
      } else {
        
        notInOld <-  
          filter(kld, year == 1990 + 23 & ticker == bad23$ticker[i]) %>% 
          distinct(name) %>% 
          pull(name)
        
      }
      
      # create new ticker names for the duplicated new firms
      set.seed(61802 + 23 + i)
      newTicks <- 
        foreach (j = 1:length(notInOld), .combine = c) %do%
        {
          
          paste(sample(letters, size = 4, replace = TRUE), collapse = "")
          
        }
      rm(j)
      newTicks <- paste(bad23$ticker[i], newTicks, sep = "~~")
      
      # take note of the new tickers for future use
      newTicksLogI <- 
        tibble(
          
          name         = notInOld,
          oldTick      = bad23$ticker[i],
          newTick      = newTicks,
          yearCreated  = 1990 + 23,
          yearLastUsed = 1990 + 23
          
        )
      
      if (i == 1) newTicksLog <- NULL
      newTicksLog <- bind_rows(newTicksLog, newTicksLogI)
      
      # execute the replacement
      for (j in 1:nrow(newTicksLogI))
      {
        
        kld$ticker[kld$year == (1990 + 23) & kld$name == newTicksLogI$name[j]] <- 
          newTicksLogI$newTick[j]
        
      }
      rm(j)

      
    }
    rm(bad23, newTicksLogI, oldI, i, newNamesI, newTicks, notInOld)
    firms_in_year[[23]] <- 
      filter(kld, year == 1990 + 23) %>% 
      pull(ticker) %>% 
      sort(.)
    # fixed
    
    # 24 is a bad year
    filter(kld, year == 1990 + 24) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num)) -> bad24
    # same idea
    # roll new tickers through
    for (i in 1:nrow(newTicksLog))
    {
      
      check <- filter(kld, year == 1990 + 24 & name == newTicksLog$name[i])
      
      if (nrow(check) > 0)
      {
        
        kld$ticker[kld$year == 1990 + 24 & kld$name == newTicksLog$name[i]] <- 
          newTicksLog$newTick[i]
        
        newTicksLog$yearLastUsed <- 1990 + 24
        
      }
      
      rm(check)

    }
    rm(i)
    # check again 
    filter(kld, year == 1990 + 24) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num)) -> bad24
    # better, but not good yet
    for (i in 1:nrow(bad24))
    {
      
      # see if the ticker has been represented earlier.
      oldI <- filter(kld, year < 1990 + 24 & ticker == bad24$ticker[i]) %>%
        select(name, ticker, year)
      
      # figure out which of the duplicated tickers are new
      if (nrow(oldI) > 0)
      {
        
        newNamesI <- 
          filter(kld, year == 1990 + 24 & ticker == bad24$ticker[i]) %>% 
          distinct(name) %>% 
          pull(name)
        
        notInOld <- setdiff(newNamesI, oldI$name)
        
      } else {
        
        notInOld <-  
          filter(kld, year == 1990 + 24 & ticker == bad24$ticker[i]) %>% 
          distinct(name) %>% 
          pull(name)
        
      }
      
      # create new ticker names for the duplicated new firms
      set.seed(61802 + 24 + i)
      newTicks <- 
        foreach (j = 1:length(notInOld), .combine = c) %do%
        {
          
          paste(sample(letters, size = 4, replace = TRUE), collapse = "")
          
        }
      rm(j)
      newTicks <- paste(bad24$ticker[i], newTicks, sep = "~~")
      
      # take note of the new tickers for future use
      newTicksLogI <- 
        tibble(
          
          name         = notInOld,
          oldTick      = bad24$ticker[i],
          newTick      = newTicks,
          yearCreated  = 1990 + 24,
          yearLastUsed = 1990 + 24
          
        )
      
      newTicksLog <- bind_rows(newTicksLog, newTicksLogI)
      
      # execute the replacement
      for (j in 1:nrow(newTicksLogI))
      {
        
        kld$ticker[kld$year == (1990 + 24) & kld$name == newTicksLogI$name[j]] <- 
          newTicksLogI$newTick[j]
        
      }
      rm(j)

    }
    rm(bad24, newTicksLogI, oldI, i, newNamesI, newTicks, notInOld)
    # all that's left is a double-counted one like in the earlier bad years.
    firms_in_year[[24]] <- 
      filter(kld, year == 1990 + 24) %>% 
      pull(ticker) %>% 
      sort(.)
    firms_in_year[[24]] <- sort(unique(firms_in_year[[24]]))
    
    # 25 is a bad year
    filter(kld, year == 1990 + 25) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num))
    # for some reason Homeaway is counted as Expedia.
    kld$ticker[kld$name == "HOMEAWAY, INC."] <- "AWAY"
    firms_in_year[[25]] <- sort(unique(firms_in_year[[25]]))
    
    # 27 is a bad year
    filter(kld, year == 1990 + 27) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num))
    # another extensive one.
    # roll new tickers throguh
    for (i in 1:nrow(newTicksLog))
    {
      
      check <- filter(kld, year == 1990 + 27 & name == newTicksLog$name[i])
      
      if (nrow(check) > 0)
      {
        
        kld$ticker[kld$year == 1990 + 27 & kld$name == newTicksLog$name[i]] <- 
          newTicksLog$newTick[i]
        
        newTicksLog$yearLastUsed <- 1990 + 27
        
      }
      
      rm(check)
      
    }
    rm(i)
    # check again 
    filter(kld, year == 1990 + 27) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num)) -> bad27
    # fix what's left
    for (i in 1:nrow(bad27))
    {
      
      # see if the ticker has been represented earlier.
      oldI <- filter(kld, year < 1990 + 27 & ticker == bad27$ticker[i]) %>%
        select(name, ticker, year)
      
      # figure out which of the duplicated tickers are new
      if (nrow(oldI) > 0)
      {
        
        newNamesI <- 
          filter(kld, year == 1990 + 27 & ticker == bad27$ticker[i]) %>% 
          distinct(name) %>% 
          pull(name)
        
        notInOld <- setdiff(newNamesI, oldI$name)
        
      } else {
        
        notInOld <-  
          filter(kld, year == 1990 + 27 & ticker == bad27$ticker[i]) %>% 
          distinct(name) %>% 
          pull(name)
        
      }
      
      # create new ticker names for the duplicated new firms
      set.seed(61802 + 27 + i)
      newTicks <- 
        foreach (j = 1:length(notInOld), .combine = c) %do%
        {
          
          paste(sample(letters, size = 4, replace = TRUE), collapse = "")
          
        }
      rm(j)
      newTicks <- paste(bad27$ticker[i], newTicks, sep = "~~")
      
      # take note of the new tickers for future use
      newTicksLogI <- 
        tibble(
          
          name         = notInOld,
          oldTick      = bad27$ticker[i],
          newTick      = newTicks,
          yearCreated  = 1990 + 27,
          yearLastUsed = 1990 + 27
          
        )
      
      newTicksLog <- bind_rows(newTicksLog, newTicksLogI)
      
      # execute the replacement
      for (j in 1:nrow(newTicksLogI))
      {
        
        kld$ticker[kld$year == (1990 + 27) & kld$name == newTicksLogI$name[j]] <- 
          newTicksLogI$newTick[j]
        
      }
      rm(j)
      
    }
    rm(bad27, newTicksLogI, oldI, i, newNamesI, newTicks, notInOld)
    # fixed
    firms_in_year[[27]] <- 
      filter(kld, year == 1990 + 27) %>% 
      pull(ticker) %>% 
      sort(.)
    
    # 28 is a bad year
    filter(kld, year == 1990 + 28) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num))
    # another extensive one.
    # roll new tickers throguh
    for (i in 1:nrow(newTicksLog))
    {
      
      check <- filter(kld, year == 1990 + 28 & name == newTicksLog$name[i])
      
      if (nrow(check) > 0)
      {
        
        kld$ticker[kld$year == 1990 + 28 & kld$name == newTicksLog$name[i]] <- 
          newTicksLog$newTick[i]
        
        newTicksLog$yearLastUsed <- 1990 + 28
        
      }
      
      rm(check)
      
    }
    rm(i)
    filter(kld, year == 1990 + 28) %>%
      group_by(ticker) %>% 
      summarise(num = n()) %>% 
      filter(num > 1) %>% 
      arrange(desc(num)) -> bad28
    # fix what's left
    for (i in 1:nrow(bad28))
    {
      
      # see if the ticker has been represented earlier.
      oldI <- filter(kld, year < 1990 + 28 & ticker == bad28$ticker[i]) %>%
        select(name, ticker, year)
      
      # figure out which of the duplicated tickers are new
      if (nrow(oldI) > 0)
      {
        
        newNamesI <- 
          filter(kld, year == 1990 + 28 & ticker == bad28$ticker[i]) %>% 
          distinct(name) %>% 
          pull(name)
        
        notInOld <- setdiff(newNamesI, oldI$name)
        
      } else {
        
        notInOld <-  
          filter(kld, year == 1990 + 28 & ticker == bad28$ticker[i]) %>% 
          distinct(name) %>% 
          pull(name)
        
      }
      
      # create new ticker names for the duplicated new firms
      set.seed(61802 + 28 + i)
      newTicks <- 
        foreach (j = 1:length(notInOld), .combine = c) %do%
        {
          
          paste(sample(letters, size = 4, replace = TRUE), collapse = "")
          
        }
      rm(j)
      newTicks <- paste(bad28$ticker[i], newTicks, sep = "~~")
      
      # take note of the new tickers for future use
      newTicksLogI <- 
        tibble(
          
          name         = notInOld,
          oldTick      = bad28$ticker[i],
          newTick      = newTicks,
          yearCreated  = 1990 + 28,
          yearLastUsed = 1990 + 28
          
        )
      
      newTicksLog <- bind_rows(newTicksLog, newTicksLogI)
      
      # execute the replacement
      for (j in 1:nrow(newTicksLogI))
      {
        
        kld$ticker[kld$year == (1990 + 28) & kld$name == newTicksLogI$name[j]] <- 
          newTicksLogI$newTick[j]
        
      }
      rm(j)
      
    }
    rm(bad28, newTicksLogI, oldI, i, newNamesI, newTicks, notInOld)
    firms_in_year[[28]] <- 
      filter(kld, year == 1990 + 28) %>% 
      pull(ticker) %>% 
      sort(.)
    
    # check to see if all tickers are distinct within year
    lapply(firms_in_year, all_distinct) %>% 
      unlist() %>% 
      all(.) # we're good
    
    # get rid of duplicated rows
    kld <- 
      kld %>% 
      unite(id, ticker, year, remove = FALSE) %>% 
      distinct(id, .keep_all = TRUE) %>% 
      select(-id)
    
    # useful to have this data for merging
    write_csv(x    = kld,
              file = "firm-year-level.csv")

#   ____________________________________________________________________________
#   execute the reshape                                                     ####
    
    allTicks <- sort(unique(kld$ticker))
    
    for (i in 1:length(unique(kld$year)))
    {
      
      # get the data we have.
      dati_in <- 
        kld %>% 
        filter(year == 1990 + i & ticker %in% firms_in_year[[i]]) %>% 
        select(ticker, mets_by_year[[i]]) %>% 
        rename_with(.fn   = ~ paste(., 1990 + i, sep = "_"),
                    .cols = -ticker) %>% 
        arrange(ticker)
      
      dati_out <- 
        tibble(
          
          ticker = unique(sort(kld$ticker[!(kld$ticker %in% 
                                              firms_in_year[[i]])]))
          
        )
      
      dati <- bind_rows(dati_in, dati_out) %>% 
        arrange(ticker) %>% 
        select(-ticker)
      yrsi <- rep(i, ncol(dati))
      rm(dati_in, dati_out)
      
      if (i == 1) dat <- dati
      if (i > 1)  dat <- bind_cols(dat, dati)
      
      if (i == 1) yrs <- yrsi
      if (i > 1)  yrs <- c(yrs, yrsi)
      
      rm(yrsi, dati)

    }
    
#   ____________________________________________________________________________
#   triple check it's 0/1/NA                                                ####
    
    trigger <- function (x)
    {
      
      any(!(is.na(x) | x %in% c(0,1)))
      
    }

    dat <- 
      dat %>% 
      mutate(across(.cols = everything(),
                    .fns  = as.integer))
    
    bad <- summarise(.data = dat,
              across(.cols = everything(),
                     .fns  = trigger)) %>% 
      pivot_longer(cols      = everything(),
                   names_to  = "metric",
                   values_to = "trigger") %>% 
      filter(trigger) %>% 
      pull(metric)
    
    dat$PRO_con_A_1995[dat$PRO_con_A_1995 > 1] <- 1
    dat$NUC_con_A_1999[dat$NUC_con_A_1999 > 1] <- 1

#   ____________________________________________________________________________
#   save for later                                                          ####
    
    rm(firms_in_year, kld, mets_by_year, newTicksLog, i, bad)
    save.image(file = "data-reshaped.RData")

#   end ## (not run)