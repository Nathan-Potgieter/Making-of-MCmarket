load("data/SNP_data.Rda")
source("code/impute_missing_returns.R")

library(pacman)
p_load(tidyverse, furrr, PerformanceAnalytics, tbl2xts, rmsfuns, lubridate)

# Imputing missing values
SNP_returns <- left_join(
    SNP_data %>%
        group_by(symbol) %>%
        arrange(date) %>%
        mutate(return = log(adjusted/dplyr::lag(adjusted))) %>%
        dplyr::filter(date>first(date)) %>%
        select(date, symbol, return) %>%
        spread(symbol, return) %>%
        impute_missing_returns(impute_returns_method = "Drawn_Distribution_Own") %>% #Imputing missing returns
        gather(symbol, return, -date) %>%
        ungroup(),
    SNP_data %>%
        select(date, symbol, volume, adjusted),
    by = c("date", "symbol")
)
#Freeing up memory
rm(SNP_data)
gc()

#This function generates random portfolios with option to supply sharp ratio
gen_random_port <- function(dim = 100, sharp = TRUE){
    #list of Assets from which sample
    symbols <- SNP_returns %>%
        dplyr::filter(date==first(date)) %>%
        pull(symbol)
    dim <- dim
    sample_symbols <- sample(symbols, dim)

    #Dates from which to sample
    dates <- SNP_returns %>%
        dplyr::filter(symbol == "A") %>%
        select(date)
    indx <- sample.int(nrow(dates) - 251, 1)
    start_date <- dates[indx,]
    end_date <- dates[indx + 251,]
    sample_dates <- dates %>% dplyr::filter(date >= start_date[[1]] &
                   date <= end_date[[1]]) %>% pull()

    if(sharp == FALSE){
        training_data <- list(dates = list(sample_dates),
                              symbols = list(sample_symbols)) %>% as_tibble()
    }else
        if(sharp == TRUE){
            #setting rebalance months
            RebMonths <- c(1,4,7,10)
            EQweights <-
                SNP_returns %>%
                dplyr::filter(symbol %in% sample_symbols &
                                  date %in% sample_dates) %>%
                select(date, symbol, return) %>%
                mutate(Months = as.numeric(format(date, format = "%m")),
                       YearMonths = as.numeric(format(date, format = "%Y%m"))) %>%
                dplyr::filter(Months %in% RebMonths) %>%
                group_by(YearMonths, Months, symbol) %>%
                dplyr::filter(date == last(date)) %>%
                ungroup() %>%
                group_by(date) %>%
                mutate(weight = 1/n()) %>%
                select(date, symbol, weight) %>%
                spread(symbol, weight) %>% ungroup() %>%
                tbl_xts()

            # Return wide data
            Returns <-
                SNP_returns %>%
                dplyr::filter(symbol %in% sample_symbols &
                                  date %in% sample_dates) %>%
                select(date, symbol, return) %>%
                spread(symbol, return) %>% tbl_xts()

            #calculating portfolio returns
            EW_RetPort <- rmsfuns::Safe_Return.portfolio(Returns,
                                               weights = EQweights, lag_weights = TRUE,
                                               verbose = TRUE, contribution = TRUE,
                                               value = 100, geometric = TRUE)
            Sharp <- SharpeRatio(EW_RetPort$returns, FUN="StdDev", annualize = TRUE)

            training_data <- list(dates = list(sample_dates),
                                  symbols = list(sample_symbols),
                                  sharp = Sharp[1]) %>% as_tibble()
    }

}
#------------------------------------------
#First I generate the labeled training data
#------------------------------------------

#Generating N random portfolios, with Sharp ratio's.
plan(multiprocess)
training_data_sharp <-
    1:100 %>% future_map(~gen_random_port(dim = 50), .progress = TRUE) %>% reduce(bind_rows)

#saving and loading created data.
save(training_data_sharp, file = "data/training_data_sharp.Rda")
load("data/training_data_sharp.Rda")

# Separating by sharp ratio into market types.
stressed_market <- training_data_sharp %>% dplyr::filter(sharp < -0.05) %>% select(-sharp)

rally_market <- training_data_sharp %>% dplyr::filter(sharp > 0.2) %>% select(-sharp)

normal_market <- training_data_sharp %>% dplyr::filter(sharp >= -0.5 & sharp <= 0.2) %>% select(-sharp)

#This function gathers the SNP data corresponding to the portfolios described above
get_market_data <- function(index_df, i){
    SNP_returns %>%
        select(date, symbol, return) %>%
            dplyr::filter(symbol %in% index_df$symbols[[i]] &
                       date %in% index_df$dates[[i]])
}

#Separating portfolio by sharp ratio
stressed_market_data <- 1:nrow(stressed_market) %>% map(~get_market_data(stressed_market, .x))

rally_market_data <- 1:nrow(rally_market) %>% map(~get_market_data(rally_market, .x))

normal_market_data <- 1:nrow(normal_market) %>% map(~get_market_data(normal_market, .x))

#Calculating correlations
calc_cor <- function(df, i, method = "spearman"){
    df[[i]] %>% select(date, symbol, return) %>%
        group_by(symbol) %>%
        mutate(return = (return - mean(return, na.rm = TRUE))/sd(return, na.rm = TRUE)) %>% #normalizing returns
        spread(symbol, return) %>% select(-date) %>% cor(method = method)
}


stressed_market_corr <- 1:length(stressed_market_data) %>% map(~calc_cor(stressed_market_data, .x))

rally_market_corr <- 1:length(rally_market_data) %>% map(~calc_cor(rally_market_data, .x))

normal_market_corr <- 1:length(normal_market_data) %>% map(~calc_cor(normal_market_data, .x))

#Joining and saving datasets
labeled_training_data <- list(stressed_market = stressed_market_corr,
                              rally_market = rally_market_corr,
                              normal_market = normal_market_corr)
save(labeled_training_data, file = "data/labeled_training_data.Rda")

rm()
gc()

#------------------------------------------
#Now to generate unlabeled training data
#------------------------------------------

plan(multiprocess)
training_data_indx <-
    1:10000 %>% future_map(~gen_random_port(dim = 50, sharp = FALSE), .progress = TRUE) %>% reduce(bind_rows)

market_data <- 1:nrow(training_data) %>% future_map(~get_market_data(training_data, .x), .progress = TRUE)

training_data <- 1:length(market_data) %>% map(~calc_cor(market_data, .x))

#saving dataframes
save(market_data, file = "data/market_data.Rda")
save(training_data, file = "data/training_data.Rda")
save(training_data, file = "data/training_data.csv")
