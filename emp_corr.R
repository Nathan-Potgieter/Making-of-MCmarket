pacman::p_load(tidyverse, tidyquant)
load("data/SNP_data.Rda")

SNP_data

date <- Sys.Date()
Tickers <- tq_index("SP500") %>%  slice_max(order_by = weight, n = 100) %>% pull(symbol)

get_data <- function(ticker = "AAPL", from) {
    df <- tq_get(ticker, from = from) %>% mutate(symbol = rep(ticker, length(date)))
    return(df)
}

SNP_top100_data <-
    1:length(stock_list) %>% map(~get_data(ticker = Tickers[[.]], from = date)) %>% bind_rows()
save(SNP_top100_data, file = "data/SNP_top100_data.rda")

get_data(ticker = "AAPL", from = as.Date("2015-01-01"))

tidyquant::tq_get_options()