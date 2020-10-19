# Note I set impute_returns_method to NONE by default. This will produce a warning if there is indeed a missing return in the matrix...

impute_missing_returns <- function(return_mat, impute_returns_method = "NONE", Seed = NULL){
    # Make sure we have a date column called date:
    if( !"date" %in% colnames(return_mat) ) stop("No 'date' column provided in return_mat. Try again please.")

    # Note my use of 'any' below...
    # Also note that I 'return' return_mat - which stops the function and returns return_mat.
    if( impute_returns_method %in% c("NONE", "None", "none") ) {
        if( any(is.na(return_mat)) ) warning("There are missing values in the return matrix.. Consider maybe using impute_returns_method = 'Drawn_Distribution_Own' / 'Drawn_Distribution_Collective'")
        return(return_mat)
    }


    if( impute_returns_method  == "Average") {

        return_mat <-
            return_mat %>% gather(Stocks, Returns, -date) %>%
            group_by(date) %>%
            mutate(Avg = mean(Returns, na.rm=T)) %>%
            mutate(Avg = coalesce(Avg, 0)) %>% # date with no returns - set avg to zero
            ungroup() %>%
            mutate(Returns = coalesce(Returns, Avg)) %>% select(-Avg) %>% spread(Stocks, Returns)

        # That is just so much easier when tidy right? See how I gathered and spread again to give back a wide df?

    } else

        if( impute_returns_method  == "Drawn_Distribution_Own") {

            set.seed(Seed)
            N <- nrow(return_mat)
            return_mat <-

                left_join(return_mat %>% gather(Stocks, Returns, -date),
                          return_mat %>% gather(Stocks, Returns, -date) %>% group_by(Stocks) %>%
                              do(Dens = density(.$Returns, na.rm=T)) %>%
                              ungroup() %>% group_by(Stocks) %>% # done to avoid warning.
                              do(Random_Draws = sample(.$Dens[[1]]$x, N, replace = TRUE, prob=.$Dens[[1]]$y)),
                          by = "Stocks"
                ) %>%  group_by(Stocks) %>% mutate(Row = row_number()) %>% mutate(Returns = coalesce(Returns, Random_Draws[[1]][Row])) %>%
                select(-Random_Draws, -Row) %>% ungroup() %>% spread(Stocks, Returns)

        } else

            if( impute_returns_method  == "Drawn_Distribution_Collective") {

                set.seed(Seed)
                NAll <- nrow(return_mat %>% gather(Stocks, Returns, -date))

                return_mat <-
                    bind_cols(
                        return_mat %>% gather(Stocks, Returns, -date),
                        return_mat %>% gather(Stocks, Returns, -date) %>%
                            do(Dens = density(.$Returns, na.rm=T)) %>%
                            do(Random_Draws = sample(.$Dens[[1]]$x, NAll, replace = TRUE, prob=.$Dens[[1]]$y)) %>% unnest(Random_Draws)
                    ) %>%
                    mutate(Returns = coalesce(Returns, Random_Draws)) %>% select(-Random_Draws) %>% spread(Stocks, Returns)

            } else

                if( impute_returns_method  == "Zero") {
                    warning("This is probably not the best idea but who am I to judge....")
                    return_mat[is.na(return_mat)] <- 0

                } else
                    stop("Please provide a valid impute_returns_method method. Options include:\n'Average', 'Drawn_Distribution_Own', 'Drawn_Distribution_Collective' and 'Zero'.")
    return(return_mat)
}
