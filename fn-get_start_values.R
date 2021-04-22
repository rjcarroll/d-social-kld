#
#
#    name: fn-get_start_values.R
#      by: rjc
#    date: 2021-04-23
#    what: this gets starting values
#
#

#   ____________________________________________________________________________
#   there's just the one                                                    ####

    # a function that makes something like a summated rating scale.
    get_start_values <-
      
      function(dat)
      {

        # convert concerns to minus 1 instead of 1
        dat <- 
          dat %>%
          mutate(across(.cols = contains("_con_"),
                        .fns  = ~ . * -1))
        
        # get NAs per row
        nas_per_row <- apply(dat, 1, function(x) sum(is.na(x)))
        
        # get NAs per column
        nas_per_column <- apply(dat, 2, function(x) sum(is.na(x)))
        
        # make the output
        res <- vector(mode = "list")
        res$firm <- apply(dat, 1, sum, na.rm = T)
        res$firmper <- res$firm/nas_per_row
        res$met <- apply(dat, 2, sum, na.rm = T)
        res$metper <- res$met/nas_per_column
        
        res

      }

#   end ## not run