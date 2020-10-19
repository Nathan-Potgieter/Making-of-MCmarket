load("data/SNP_returns.Rda")
#note that this first step takes a few minuets.
garch_high_vol <-
    SNP_returns %>%
    select(date, symbol, return) %>%
    dplyr::filter(symbol %in% high_vol) %>%
    split(.$symbol) %>% length()
    map(~garchFit(formula = ~arma(1)+aparch(1,1), cond.dist = "sstd", data = .$return, trace = FALSE))

garch_low_vol <-
    SNP_returns %>%
    select(date, symbol, return) %>%
    dplyr::filter(symbol %in% low_vol) %>%
    split(.$symbol) %>%
    map(~garchFit(formula = ~arma(1)+aparch(1,1), cond.dist = "sstd", data = .$return, trace = FALSE))

garch_medium_vol <-
    SNP_returns %>%
    select(date, symbol, return) %>%
    dplyr::filter(symbol %in% medium_vol) %>%
    split(.$symbol) %>%
    map(~garchFit(formula = ~arma(1)+aparch(1,1), cond.dist = "sstd", data = .$return, trace = FALSE))


#Creating a dataframe containing the GARCH residues.
get_resid <- function(garch) {

    residuals <- suppressMessages(
        1:length(garch) %>% map(~residuals(garch[[.]])) %>%
            reduce(bind_cols)
    )
    colnames(residuals) <- names(garch)
    return(residuals)
}

resid_high_vol <- get_resid(garch_high_vol)
resid_low_vol <- get_resid(garch_low_vol)
resid_medium_vol <- get_resid(garch_medium_vol)

#Creating a dataframe containing GARCH parameter estimates.
get_coef <- function(garch) {

    coef <- 1:length(garch) %>% map(~coef(garch[[.]])) %>%
        reduce(bind_rows) %>%
        mutate(Ticker = names(garch), .before = mu)

    return(coef)
}

coef_high_vol <- get_coef(garch_high_vol)
coef_low_vol <- get_coef(garch_low_vol)
coef_medium_vol <- get_coef(garch_medium_vol)

save(list = c("resid_high_vol",
              "resid_low_vol",
              "resid_medium_vol",
              "coef_high_vol",
              "coef_low_vol",
              "coef_medium_vol"),
     file = "data/garch.Rda")
#rm(garch)    #unquote later too save memory



