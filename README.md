README
================

# Simulating Asset Returns

This is the README for Nathan Potgieter’s financial econometrics
project, in which a framework the Monte Carlo simulation of asset
markets is developed.

## Aim

The aim of this project is to develop a general and easy to use Monte
Carlo simulation package that generates asset return data, with a
prespecified correlation structure and dynamic dependencies. Ideally the
user will be able to adjust a “leverage” parameter, which will determine
the markets left-tail dependency, and in turn effect likelihood of
entering a “crisis period” characterized by extreme joint drawdowns.

Elliptical copulas are used to induce the correlation in the simulated
data, while Archmedian copulas are used induce greater left-tail
dependencies.

The data will also be simulated to exhibit volatility clustering, this
is accomplished by utilizing an ARIMA(p,q) + APGARCH(q,p) model, the
parameters of which can be adjusted to induce alternative risk
characteristics. Various ARIMA(p,q) + APGARCH(q,p) structures can be
called on to induce mean and variance persistence.

## Monte Carlo Framework

The Monte Carlo simulation routine involves the following steps:

This example generates 252 days worth of returns, for 50 Assets across N
markets.

1.  Draw a series of 252 random uniformly distributed numbers
    (corresponding to 252 trading days), across a set of 50 variables
    (or 50 assets), from a multivariate distribution with a given
    correlation matrix.
    -   This is accomplished using Euclidean (Gaussian or t-copula) and
        Archmediean copula’s and can easily be done using the rcopula
        function.
2.  Convert the uniformly distributed marginal distributions into
    something that more resembles the distribution of asset returns. For
    example one could convert them into normal or skewed-generalized t
    distributions.
    -   This is done the same way one would convert p-values into test
        statistics using the dnorm() and dsgt() functions respectively.
3.  The next step is to induce mean and variance persistence to the
    series, by plugging them into a ARIMA(p,q) + GARCH(q,p) equation as
    innovations.
    -   If the parameters are set accordingly the resulting series
        should resemble asset returns.
4.  The final step is to repeat the first 3 steps N times to generate an
    ensemble of asset markets, each with the same prespecified
    structure.

#### Loading Packages

``` r
library(pacman)
p_load(tidyverse, copula, fGarch, lubridate, forecast, bizdays, sgt, glue)
p_load(tbl2xts)
```

# The Set up

## Generating Ad Hoc Covarience matrix

In this section I developed a simple function that allows the user to
easily generate a correlation matrix with a desired cluster structure.
This will be used as a key input when simulating our financial markets.
Note that the majority of the code was written by Nico Katzke. The
function is located in the gen\_corr.R code file.

### gen\_corr’s arguments

1.  N - is the number of assets in the universe

2.  Clusters - a character string specifying the type of cluster
    structure. Available options are “none”, for a correlation matrix
    with no clusters, “non-overlapping” for a correlation matrix with
    number one layer of clusters, and “overlapping” for a correlation
    matrix with Num\_Layers and Num\_clusters per layer.

3.  Num\_Clusters - if Clusters is equal to “non-overlapping” or “none”
    then Num\_Clusters is an integer value specifying the number of
    clusters. If Clusters = “overlapping” then Num\_Clusters must be a
    vector of length equal to Num\_Layers specifying the number of
    clusters per layer.

4.  Num\_Layers - an integer value between 1 and 4, specifying the
    number of cluster layers. Only needed of using “overlapping”
    clusters.

``` r
#Co-Varience matrix generatimg function

gen_corr <- 
  function(N = 50, Clusters = c("none", "non-overlapping", "overlapping") , Num_Clusters = NULL, Num_Layers = NULL){
    
Grps <- Num_Clusters
#set.seed(123)
    
if(Clusters == "none"){
    # Unclustered covariance matrix
    Sigma <- diag(N)
    for (i in 1:N) for (j in 1:N) Sigma[i,j] <- 0.9^abs(i-j)
    Sigma <- propagate::cor2cov(Sigma, runif(N, 1, 5))
    corr <- cov2cor(Sigma)
} else

if(Clusters == "non-overlapping"){
    #----------------------
    # distinct non-overlapping clusters:
    #----------------------
    
    if(is.null(Num_Clusters)) stop("Please provide a valid Num_Clusters argument when using Overlapping clusters")
    
    
    Sigma <- matrix(0.9, N, N)
    diag(Sigma) <- 1

    
for (i in 1:Grps) {
      ix <- seq((i-1) * N / Grps + 1, i * N / Grps)
      Sigma[ix, -ix] <- 0.0001                       #think about
    }
    Sigma <- propagate::cor2cov(Sigma, runif(N, 1, 5))
    corr <- cov2cor(Sigma)
} else
  
if(Clusters == "overlapping"){
    #----------------------
    # distinct overlapping clusters:
    #----------------------
  
  if(is.null(Num_Layers)|Num_Layers<2){
      stop("Please provide a valid Num_Layers argument when using Overlapping clusters")
      }else
  if(length(Num_Clusters) != Num_Layers){
      stop("Please provide a Num_Clusters argument with length equal to Num_Layers")
  }
    
  
    Sigma <- matrix(0.9, N, N)
    diag(Sigma) <- 1

    for (i in 1:Grps[1]) {
      ix <- seq((i-1) * N / Grps[1] + 1, i * N / Grps[1])
      Sigma[ix, -ix] <- 0.7
    }
    if(Num_Layers>=2){
        for (i in 1:Grps[2]) {
          ix <- seq((i-1) * N / Grps[2] + 1, i * N / Grps[2])
          Sigma[ix, -ix] <- 0.5
        } }
    if(Num_Layers>=3){
        for (i in 1:Grps[3]) {
      ix <- seq((i-1) * N / Grps[3] + 1, i * N / Grps[3])
      Sigma[ix, -ix] <- 0.3
        } }
    if(Num_Layers>=4){
        for (i in 1:Grps[4]) {
      ix <- seq((i-1) * N / Grps[4] + 1, i * N / Grps[4])
      Sigma[ix, -ix] <- 0.15
        } } 
    }

    Sigma <- propagate::cor2cov(Sigma, runif(N, 1, 5))  #Is this necessary???
    corr <- cov2cor(Sigma)

return(corr)

  }
```

Demonstrating the use of gen\_corr

``` r
source("code/gen_corr.R")
gen_corr(N = 60, Clusters = "overlapping", Num_Layers = 4, Num_Clusters = c(10,5,3,2)) %>% ggcorrplot::ggcorrplot(title = "Overlapping Clusters", hc.order = TRUE)
```

<img src="README_files/figure-gfm/using gen_corr-1.png" width="80%" height="80%" />

## Generating a Dataset of Emperical Correlation Matrix’s

I now use S&P500 data since 1/01/2000 to sample correlation matrices
that will be used to train CorrGAN and can be used as inputs in the
simulation.

#### Getting SNP 500 data since 2000

-   See “code/SNP\_data.R” to see how the SNP\_data.Rda file was
    created.

``` r
p_load(tidyverse, tidyquant, lubridate)
source("code/impute_missing_returns.R")

# save current system date to a variable
today <- Sys.Date()
# subtract 90 months from the current date
#date <- today %m+% months(-90)
date <- as.Date("2000-01-01")

# Getting all tickers in SP500
stock_list <- tq_index("SP500") %>%
    arrange(symbol) %>%
    mutate(symbol = case_when(symbol == "BRK.B" ~ "BRK-B",
                              symbol == "BF.B" ~ "BF-B",
                              TRUE ~ as.character(symbol))) %>%
    pull(symbol)

# This function gets the data for a ticker and date
get_data <- function(ticker = "AAPL", from){
    df <- tq_get(ticker, from = from) %>% mutate(symbol = rep(ticker, length(date)))
    return(df)
}

SNP_data <-
    1:length(stock_list) %>% map(~get_data(ticker = stock_list[[.]], from = date)) %>% bind_rows()
save(SNP_data, file = "data/SNP_data.Rda")

#Calculating and imputing missing returns
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
        group_by(symbol) %>%
        mutate(sd = sd(return, na.rm = T)*sqrt(252)) %>%  # calculate annualized SD
        ungroup(),
    SNP_data %>%
        select(date, symbol, volume, adjusted),
    by = c("date", "symbol")
)
save(SNP_returns, file = "data/SNP_returns.Rda")
```

#### Partisioning SNP Data into Rally, Normal and Stressed Markets

-   See “code/get\_training\_data.R” to see how the training data sets
    were generated. The data set is build to contain 3 classes of
    correlation matrices, defined as follows:

-   ‘stressed market’: A market is ‘stressed’ whenever the equi-weighted
    basket of stocks has a Sharpe below -0.5 over the year of study (252
    trading days).

-   ‘rally market’: A market is ‘rallying’ whenever the equi-weighted
    basket of stocks under has a Sharpe above 2 over the year of study
    (252 trading days).

-   ‘normal market’: A market is ‘normal’ whenever the equi-weighted
    basket of stocks under has a Sharpe in-between -0.5 and 2 over the
    year of study (252 trading days).

Note that this methodology is consistent with that used in
<https://marti.ai/qfin/2020/02/03/sp500-sharpe-vs-corrmats.html>.

``` r
load("data/SNP_data.Rda")
source("code/impute_missing_returns.R")

library(pacman)
p_load(tidyverse, furrr, PerformanceAnalytics, tbl2xts, rmsfuns, lubridate, fitHeavyTail)

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
    indx <- sample.int(nrow(dates) - 252, 1)
    start_date <- dates[indx,]
    end_date <- dates[indx + 252,]
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
    1:10000 %>% future_map(~gen_random_port(dim = 50), .progress = TRUE) %>% reduce(bind_rows)
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
calc_cor <- function(df, i){
    df[[i]] %>% select(date, symbol, return) %>%
        spread(symbol, return) %>% select(-date) %>%
        fitHeavyTail::fit_mvt() %>% .$cov %>% cov2cor()
}

stressed_market_corr <- 1:length(stressed_market_data) %>% map(~calc_cor(stressed_market_data, .x))

rally_market_corr <- 1:length(rally_market_data) %>% map(~calc_cor(rally_market_data, .x))

normal_market_corr <- 1:length(normal_market_data) %>% map(~calc_cor(normal_market_data, .x))

#Joining and saving datasets
labeled_training_data <- list(stressed_market = stressed_market_corr,
                              rally_market = rally_market_corr,
                              normal_market = normal_market_corr)
save(labeled_training_data, file = "data/labeled_training_data.Rda")

rm(stressed_market)
gc()

#------------------------------------------
#Now to generate unlabeled training data
#------------------------------------------

plan(multiprocess)
training_data_indx <-
    1:10000 %>% future_map(~gen_random_port(dim = 50, sharp = FALSE), .progress = TRUE) %>% reduce(bind_rows)

market_data <- 1:nrow(training_data) %>% future_map(~get_market_data(training_data, .x), .progress = TRUE)

training_data <- 1:length(market_data) %>% map(~calc_cor(market_data, .x))


save(training_data, file = "data/training_data.Rda")
```

#### An Example of Each Correlation Matrix Type

<img src="README_files/figure-gfm/showing some training data-1.png" width="80%" height="80%" /><img src="README_files/figure-gfm/showing some training data-2.png" width="80%" height="80%" /><img src="README_files/figure-gfm/showing some training data-3.png" width="80%" height="80%" />

# Step 1: Draw a series of random uniformly distributed numbers across a set of variables with a specified dependence structure.

## Generating Random Draws with Numerous Copula Functions

### Elliptal copulas

Elliptal copulas such as the Gaussian and the student t copulas, allow
us to specify a correlation matrix before randomly selecting
observations from the multivariate distribution. Doing so allows one to
produce random draws of uniformly distributed variables, that contain
the correlation structure and joint distribution specified by the
copula. The chunk of code below demonstrates this functionality.

``` r
#loading copula package
pacman::p_load(copula)

#generating corr  matrix object
corr <- gen_corr(N = 50, Clusters = "overlapping", Num_Layers = 3, Num_Clusters = c(10, 5, 2))

#generating copula objects   
Ncop <- ellipCopula(family = "normal", dispstr = "un", param = P2p(corr), dim = 50)
Tcop <- ellipCopula(family = "t", dispstr = "un", param = P2p(corr), dim = 50)

#generating 252 random draws for each of the N variables
set.seed(123)
rn <- rCopula(copula = Ncop, n = 252)
rt <- rCopula(copula = Tcop, n = 252)

#Checking if the correlation structure was maintained
p_load(patchwork)
# Original corr
p1 <- ggcorrplot::ggcorrplot(corr, hc.order = TRUE) + 
  labs(title = "Input Correlation Matrix") +
  scale_x_discrete(labels = NULL) + scale_y_discrete(labels = NULL) +
  theme(legend.position = "bottom")
# corr from random draws form norm and t copula
p2 <- fitHeavyTail::fit_mvt(rn) %>% .$cov %>% cov2cor() %>% 
  ggcorrplot::ggcorrplot(hc.order = TRUE) + 
  labs(subtitle = "Normal Copula") +
  scale_x_discrete(labels = NULL) + scale_y_discrete(labels = NULL) +
  theme(legend.position = "none")

p3 <- fitHeavyTail::fit_mvt(rt) %>% .$cov %>% cov2cor() %>% 
  ggcorrplot::ggcorrplot(hc.order = TRUE) + 
  labs(subtitle = "T-Copula") +
  scale_x_discrete(labels = NULL) + scale_y_discrete(labels = NULL) +
  theme(legend.position = "none")

p_load(printr)
summary(rn)
```

|     | V1             | V2               | V3               | V4               | V5               | V6               | V7               | V8               | V9               | V10              | V11              | V12              | V13              | V14               | V15              | V16               | V17              | V18              | V19               | V20               | V21               | V22               | V23               | V24              | V25              | V26              | V27              | V28              | V29              | V30              | V31              | V32             | V33              | V34              | V35              | V36              | V37              | V38              | V39              | V40              | V41              | V42              | V43             | V44             | V45              | V46              | V47              | V48              | V49             | V50              |
|:----|:---------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:------------------|:-----------------|:------------------|:-----------------|:-----------------|:------------------|:------------------|:------------------|:------------------|:------------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:-----------------|:----------------|:----------------|:-----------------|:-----------------|:-----------------|:-----------------|:----------------|:-----------------|
|     | Min. :0.0040   | Min. :0.006049   | Min. :0.002107   | Min. :0.002224   | Min. :0.001711   | Min. :0.003406   | Min. :0.003609   | Min. :0.003013   | Min. :0.008587   | Min. :0.003547   | Min. :0.001026   | Min. :0.008247   | Min. :0.001639   | Min. :0.0004023   | Min. :0.003664   | Min. :0.0008716   | Min. :0.001115   | Min. :0.000027   | Min. :0.0004764   | Min. :0.0000399   | Min. :0.0003142   | Min. :0.0006477   | Min. :0.0007081   | Min. :0.002966   | Min. :0.001216   | Min. :0.003294   | Min. :0.007127   | Min. :0.005926   | Min. :0.005297   | Min. :0.005952   | Min. :0.007729   | Min. :0.01515   | Min. :0.005879   | Min. :0.001069   | Min. :0.006362   | Min. :0.001005   | Min. :0.003714   | Min. :0.002156   | Min. :0.001275   | Min. :0.006033   | Min. :0.003553   | Min. :0.006991   | Min. :0.01212   | Min. :0.01201   | Min. :0.005037   | Min. :0.005705   | Min. :0.009809   | Min. :0.005633   | Min. :0.00275   | Min. :0.006285   |
|     | 1st Qu.:0.2913 | 1st Qu.:0.260398 | 1st Qu.:0.268375 | 1st Qu.:0.259785 | 1st Qu.:0.274867 | 1st Qu.:0.274325 | 1st Qu.:0.294156 | 1st Qu.:0.288425 | 1st Qu.:0.308749 | 1st Qu.:0.298146 | 1st Qu.:0.268024 | 1st Qu.:0.298956 | 1st Qu.:0.306712 | 1st Qu.:0.2929735 | 1st Qu.:0.262164 | 1st Qu.:0.2512954 | 1st Qu.:0.265083 | 1st Qu.:0.269510 | 1st Qu.:0.2478857 | 1st Qu.:0.2522780 | 1st Qu.:0.2603263 | 1st Qu.:0.2400959 | 1st Qu.:0.2670890 | 1st Qu.:0.254333 | 1st Qu.:0.254212 | 1st Qu.:0.267642 | 1st Qu.:0.287304 | 1st Qu.:0.249182 | 1st Qu.:0.226799 | 1st Qu.:0.244830 | 1st Qu.:0.218558 | 1st Qu.:0.26162 | 1st Qu.:0.219807 | 1st Qu.:0.227548 | 1st Qu.:0.226846 | 1st Qu.:0.240695 | 1st Qu.:0.235379 | 1st Qu.:0.220862 | 1st Qu.:0.233688 | 1st Qu.:0.254933 | 1st Qu.:0.222121 | 1st Qu.:0.223913 | 1st Qu.:0.22789 | 1st Qu.:0.22098 | 1st Qu.:0.202687 | 1st Qu.:0.201772 | 1st Qu.:0.225565 | 1st Qu.:0.208706 | 1st Qu.:0.22236 | 1st Qu.:0.234473 |
|     | Median :0.5434 | Median :0.546300 | Median :0.547710 | Median :0.532551 | Median :0.543623 | Median :0.549016 | Median :0.543377 | Median :0.533298 | Median :0.533712 | Median :0.504722 | Median :0.526656 | Median :0.543051 | Median :0.518387 | Median :0.5515773 | Median :0.511259 | Median :0.5155747 | Median :0.522784 | Median :0.521189 | Median :0.5707847 | Median :0.5163259 | Median :0.5061045 | Median :0.5203811 | Median :0.5183899 | Median :0.505304 | Median :0.492845 | Median :0.525529 | Median :0.525535 | Median :0.505862 | Median :0.538242 | Median :0.504900 | Median :0.505584 | Median :0.48360 | Median :0.486611 | Median :0.469043 | Median :0.456277 | Median :0.477655 | Median :0.479151 | Median :0.473166 | Median :0.491240 | Median :0.491238 | Median :0.467227 | Median :0.455782 | Median :0.44742 | Median :0.44763 | Median :0.464029 | Median :0.449205 | Median :0.448168 | Median :0.462077 | Median :0.46800 | Median :0.441365 |
|     | Mean :0.5192   | Mean :0.515365   | Mean :0.521043   | Mean :0.524026   | Mean :0.528011   | Mean :0.522403   | Mean :0.527877   | Mean :0.529493   | Mean :0.523711   | Mean :0.519163   | Mean :0.524524   | Mean :0.519849   | Mean :0.514609   | Mean :0.5242924   | Mean :0.507877   | Mean :0.5042793   | Mean :0.510002   | Mean :0.507160   | Mean :0.5169573   | Mean :0.5104083   | Mean :0.5011769   | Mean :0.4973071   | Mean :0.5060419   | Mean :0.504966   | Mean :0.489123   | Mean :0.517274   | Mean :0.520168   | Mean :0.499533   | Mean :0.509322   | Mean :0.503921   | Mean :0.483099   | Mean :0.47766   | Mean :0.471798   | Mean :0.473442   | Mean :0.471736   | Mean :0.492342   | Mean :0.483843   | Mean :0.488183   | Mean :0.490003   | Mean :0.492984   | Mean :0.476356   | Mean :0.476640   | Mean :0.48604   | Mean :0.48300   | Mean :0.474777   | Mean :0.465217   | Mean :0.467509   | Mean :0.466425   | Mean :0.47318   | Mean :0.463435   |
|     | 3rd Qu.:0.7647 | 3rd Qu.:0.758883 | 3rd Qu.:0.778233 | 3rd Qu.:0.791575 | 3rd Qu.:0.791497 | 3rd Qu.:0.764421 | 3rd Qu.:0.740819 | 3rd Qu.:0.765897 | 3rd Qu.:0.745265 | 3rd Qu.:0.752340 | 3rd Qu.:0.763049 | 3rd Qu.:0.772789 | 3rd Qu.:0.745427 | 3rd Qu.:0.7709556 | 3rd Qu.:0.748171 | 3rd Qu.:0.7304215 | 3rd Qu.:0.759056 | 3rd Qu.:0.758749 | 3rd Qu.:0.7776007 | 3rd Qu.:0.7709514 | 3rd Qu.:0.7552242 | 3rd Qu.:0.7492251 | 3rd Qu.:0.7458391 | 3rd Qu.:0.754895 | 3rd Qu.:0.729601 | 3rd Qu.:0.781333 | 3rd Qu.:0.785622 | 3rd Qu.:0.771182 | 3rd Qu.:0.790040 | 3rd Qu.:0.756660 | 3rd Qu.:0.711921 | 3rd Qu.:0.70084 | 3rd Qu.:0.673059 | 3rd Qu.:0.702203 | 3rd Qu.:0.716678 | 3rd Qu.:0.751807 | 3rd Qu.:0.717056 | 3rd Qu.:0.735565 | 3rd Qu.:0.737668 | 3rd Qu.:0.740199 | 3rd Qu.:0.731425 | 3rd Qu.:0.702319 | 3rd Qu.:0.76009 | 3rd Qu.:0.72805 | 3rd Qu.:0.719454 | 3rd Qu.:0.721762 | 3rd Qu.:0.693500 | 3rd Qu.:0.681494 | 3rd Qu.:0.68487 | 3rd Qu.:0.684691 |
|     | Max. :0.9960   | Max. :0.998757   | Max. :0.996542   | Max. :0.993583   | Max. :0.992648   | Max. :0.996512   | Max. :0.996105   | Max. :0.993313   | Max. :0.997094   | Max. :0.997355   | Max. :0.995868   | Max. :0.997266   | Max. :0.996040   | Max. :0.9979606   | Max. :0.991795   | Max. :0.9913697   | Max. :0.996826   | Max. :0.990676   | Max. :0.9910773   | Max. :0.9947191   | Max. :0.9978228   | Max. :0.9993643   | Max. :0.9988181   | Max. :0.998551   | Max. :0.998534   | Max. :0.997108   | Max. :0.989301   | Max. :0.995131   | Max. :0.990259   | Max. :0.981684   | Max. :0.995668   | Max. :0.99902   | Max. :0.997481   | Max. :0.995905   | Max. :0.995168   | Max. :0.992198   | Max. :0.995119   | Max. :0.998803   | Max. :0.994537   | Max. :0.996372   | Max. :0.999231   | Max. :0.997015   | Max. :0.99876   | Max. :0.99755   | Max. :0.997869   | Max. :0.996131   | Max. :0.998389   | Max. :0.998085   | Max. :0.99592   | Max. :0.999283   |

``` r
p1
```

<img src="README_files/figure-gfm/Elliptal Copulas-1.png" width="80%" height="80%" />

``` r
#Notice that the underlying correlation structure has, for the most part, been maintained.
# Some Noise has been introduced
(p2+p3) + plot_annotation(title = "Output Correlation Matrices") +
  plot_layout(guides='collect') &
  theme(legend.position='bottom')  
```

<img src="README_files/figure-gfm/Elliptal Copulas-2.png" width="80%" height="80%" />

### Archimedean Copulas

Unfortunately, Elliptal copulas cannot be calibrated to exhibit varying
co-movements within the tails of the distribution. Therefore, in this
section we examine some properties of Archimedean copulas.

Archimedean copulas such as the Clayton, Frank, Gumbel and Joe exhibit
increased dependence at the tails of the multivariate distribution. In
this section we will examine the Clayton copula due to it exhibiting
enhanced left-tail dependencies. Other copulas will not be examined,
since the Clayton copula is currently the only Archimedean copula, in
the copula package, that allows random sampling from multivariate
distributions with Dim &gt; 2.

We will also have a look at the hybrid BB1-BB6 which in which exhibit
increased dynamic dependencies in both tails.

``` r
# first look at at dim=2 to get understanding of what parameter tning does

#Clayton Copula
claycop <- archmCopula(family = "clayton", param = 2, dim = 2)
Ncop <- ellipCopula(family = "normal", dispstr = "un", param = 0.5, dim = 2)
# rCopula(251, claycop)

#note how left tail dependence increases with the parameter
persp(claycop, main = "Clayon Copula" , dCopula, zlim = c(0, 15), theta = 18)
```

<img src="README_files/figure-gfm/Archimedean copula-1.png" width="80%" height="80%" />

``` r
persp(Ncop, main = "Normal Copula" , dCopula, zlim = c(0, 15), theta = 18)
```

<img src="README_files/figure-gfm/Archimedean copula-2.png" width="80%" height="80%" />

``` r
# Compairing Clayton with Normal copula
bind_rows(rCopula(5000, copula = claycop) %>% as_tibble() %>% mutate(copula = "claycop"),
          rCopula(5000, copula = Ncop) %>% as_tibble() %>% mutate(copula = "normal")) %>% 
  ggplot(aes(x=V1,y=V2)) +
        geom_point(alpha=0.5) +
        geom_density_2d_filled(alpha=0.7) +
        facet_wrap(~copula, nrow = 1) +
        labs(title = "2D kernal Density - Clayton vs Normal Copula",
             caption = "Clayton parameter = 2, Normal correlation = 0.5") +
        theme_bw() +
        theme(legend.position = "bottom") 
```

    ## Warning: The `x` argument of `as_tibble.matrix()` must have unique column names if `.name_repair` is omitted as of tibble 2.0.0.
    ## Using compatibility `.name_repair`.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_warnings()` to see where this warning was generated.

<img src="README_files/figure-gfm/Archimedean copula-3.png" width="80%" height="80%" />

``` r
#Note that the Gumbel and Joe copulas must be rotated 180 degrees to exhibit greater left tail dependence
#Gumbel Copula
# gumcop <- archmCopula(family = "gumbel", param = 4, dim = 2) %>% rotCopula()

#note how right tail dependence > left tail dependence; tail dependence increase with the parameter value.
# persp(gumcop, dCopula, zlim = c(0, 10))
# rCopula(1000, copula = gumcop) %>% plot()

#Joe copula
#joecop <- archmCopula(family = "joe", param = 3, dim = 2) %>% rotCopula()

#note how right tail dependence > left tail dependence;tail dependence increase with the parameter value at rate < gumbel
#persp(joecop, dCopula, zlim = c(0, 10))
#rCopula(1000, copula = joecop) %>% plot()


#looking at some hybrid copulas
#Galambos
#galcop <- evCopula(family = "galambos", param = 2, dim = 2) %>% rotCopula()
#persp(galcop, dCopula, zlim = c(0, 10))
#rCopula(1000, galcop) %>% plot()
```

## Generating Hybrid Copulas

Tawn’s (1988) Theorem: Shows that a copula is a convex set and every
convex combination of existing copula functions is again a copula.

Thus, if
![C\_1(U\_N)](https://latex.codecogs.com/png.latex?C_1%28U_N%29 "C_1(U_N)")
and
![C\_2(U\_N)](https://latex.codecogs.com/png.latex?C_2%28U_N%29 "C_2(U_N)")
are multivariate copula’s of dimension N and
![w](https://latex.codecogs.com/png.latex?w "w") is a weighting variable
between 0 and 1, then

<center>

![C(U\_N)=w.C\_1(U\_N)+(1-w)C\_2(U\_N)](https://latex.codecogs.com/png.latex?C%28U_N%29%3Dw.C_1%28U_N%29%2B%281-w%29C_2%28U_N%29 "C(U_N)=w.C_1(U_N)+(1-w)C_2(U_N)")

</center>

is a unique copula. Therefore, a hybrid copula (
![C(U\_N)](https://latex.codecogs.com/png.latex?C%28U_N%29 "C(U_N)") )
can be created by linearly weighting a Elliptical and Archimedean copula
of the same dimension.

See “Extreme Dependence Structures and the Cross-Section of Expected
Stock Returns” page 8 & 9.

## Some 2D Hybrid Copulas.

![Plot 1.](plots/plot_1.png) ![Plot 2.](plots/plot_2.png) ![Plot
3.](plots/plot_3.png) ![Plot 4.](plots/plot_4.png)

However, remember that each variable is currently uniformly distributed.
<img src="README_files/figure-gfm/unnamed-chunk-1-1.png" width="80%" height="80%" />

# Step 2: Converting the uniformly distributed variables to something that better resembles the distribution of asset returns.

## Looking at options for marginal distributions

Due to convenience, it has become standard to use the normal, or
student-t distribution when simulating asset returns.

However, after reading up on numerous possible marginal distributions, I
decided that the the Skewed generalized t distribution is the most
appropriate as it allows for the most flexibility. Note that I will also
include functionality to induce the marginals to be uniformly
distributed.

In fact, the SGT distribution nests 12 common probability distribution
functions (pdf). The tree diagram below indicates how one can set the
SGT parameters to achieve the desired pdf.

![Plot 5.](plots/SGTtree.png)

### The skewed generalizd t distribution with different parameters

The code below demonstrates how the p, q and
![\\lambda](https://latex.codecogs.com/png.latex?%5Clambda "\lambda")
functions influence the SGT distribution.
<img src="README_files/figure-gfm/unnamed-chunk-2-1.png" width="80%" height="80%" /><img src="README_files/figure-gfm/unnamed-chunk-2-2.png" width="80%" height="80%" /><img src="README_files/figure-gfm/unnamed-chunk-2-3.png" width="80%" height="80%" />

### Calibrating the SGT with Architypal Low, Medium and High Risk Assets

We now look at data on each share in the S&P500 over the last 90 months.
The shares with the top 5% highest annualized SD’s are used to model an
archetypal high risk asset, shares with the 5% lowest annualized SD’s
are used to model an archetypal low risk asset, while shares with SD
between the 45th and 55th percentile are used to model the medium risk
asset. - See “code/SNP\_data.R” to see how the SNP\_data.Rda file was
created.

``` r
load("data/SNP_returns.Rda")

high_vol <-
  SNP_returns %>% 
  dplyr::filter(date==last(date)) %>% 
  arrange(desc(sd)) %>% 
  group_by(date) %>% 
  slice_max(., order_by = sd, prop = 0.1 ) %>%  #select top 10% SD's
  pull(symbol)

low_vol <-
  SNP_returns %>% 
  dplyr::filter(date==last(date)) %>% 
  arrange(sd) %>% 
  slice_min(., order_by = sd, prop = 0.1) %>% #select bot 10% SD's
  pull(symbol)

medium_vol <-
  SNP_returns %>% 
  dplyr::filter(date==last(date)) %>% 
  arrange(sd) %>% 
  group_by(date) %>% 
  dplyr::filter(sd>=quantile(sd, probs = 0.45, na.rm = T) &
                  sd<=quantile(sd, probs = 0.55, na.rm = T)) %>% 
  pull(symbol)
```

# Plotting low, medium and high risk returns.

    ## Coordinate system already present. Adding new coordinate system, which will replace the existing one.

<img src="README_files/figure-gfm/Plotting returns-1.png" width="80%" height="80%" />

### Estimating SGT

-   See “code/filter\_resid.R” to see how filtered residuals (garch.Rda)
    was obtained.
-   See “code/estimate\_sgt.R” to see how the parameters of the sgt were
    estimated.

``` r
load("data/garch.Rda")

estimate_sgt <- function(df, start = NULL){
  x <- df[[1]]
  X.f <- X ~ x
  if(is.null(start)) start <- list(mu = 0, sigma = 0.03, lambda = -0.02, p = 1.5, q = 2.25)
  result <- sgt.mle(X.f = X.f, start = start, finalHessian = "BHHH")
  summary(result)
}

# High Vol stocks
df <- resid_high_vol %>% 
  dplyr::mutate(date = unique(SNP_returns$date)) %>% 
  gather(symbol, return, -date) %>% 
  dplyr::filter(date>first(date)) %>% 
  select(return)

#start <- list(mu = -0.0001957195, sigma = 0.04217965, lambda = -0.0062424590, p= 1.452239, q = 2.058042)
 # sgt_high_vol <- estimate_sgt(df, start = start)
 # save(sgt_high_vol, file = "data/sgt_high_vol.Rda")
load(file = "data/sgt_high_vol.Rda")

# Low Vol stocks
df <- resid_low_vol %>% 
  mutate(date = unique(SNP_returns$date)) %>% 
  gather(symbol, return, -date) %>% 
  dplyr::filter(date>first(date)) %>% 
  select(return)

 # sgt_low_vol <- estimate_sgt(df)
 # save(sgt_low_vol, file = "data/sgt_low_vol.Rda")
load(file = "data/sgt_low_vol.Rda")


#Medium Vol stocks
df <- resid_medium_vol %>% 
  mutate(date = unique(SNP_returns$date)) %>% 
  gather(symbol, return, -date) %>% 
  dplyr::filter(date>first(date)) %>% 
  select(return)

 # sgt_medium_vol <- estimate_sgt(df)
 # save(sgt_medium_vol, file = "data/sgt_medium_vol.Rda")
load(file = "data/sgt_medium_vol.Rda")
```

<img src="README_files/figure-gfm/plotting marginal distributions-1.png" width="80%" height="80%" /><img src="README_files/figure-gfm/plotting marginal distributions-2.png" width="80%" height="80%" /><img src="README_files/figure-gfm/plotting marginal distributions-3.png" width="80%" height="80%" />

## Simulating Innovations

This section introduces the sim\_inno function, which is designed to
carry out the first two steps of this Monte Carlo framework.

### sim\_inno

The function below generates randomly distributed numbers from a hybrid
t and clayton copula. Need to think about how to calibrate df and
claycop parameters.

Arguments

-   Corr this is a correlation matrix uses as the parameter for the
    elliptical copula
-   elliptal\_copula family name of elliptal copula. Default is to use
    “t”, but “norm” is also accepted
-   df\_ellip a positive integer specifying the degrees of freedom for
    the student t elliptic copula. Only required when using
    elliptal\_copula = “t”.
-   left\_cop\_param a positive integer specifying the parameter of the
    **Clayton** copula.
-   left\_cop\_weight a value between 0 and 1 corresponding to the
    weight assigned to the left copula, when generating random draws
    from a hybrid copula.
-   marginal\_dist a character string specifying the marginal
    distribution of the simulated data. Must be “norm” or “t”, with the
    default generating uniformly distributed marginals.
-   df\_marginal\_dist a positive integer specifying the degrees of
    freedom parameter of the “t” distributed marginals.

``` r
sim_inno <- function(corr, elliptal_copula = c("norm", "t"), 
                  df_ellip = NULL, left_cop_param = 5,
                  left_cop_weight = 0, T = 251, 
                  marginal_dist = NULL, marginal_dist_model = NULL, sd_md = NULL) {
  
  N <- nrow(corr)
  T <- T + 5
  Cor <- P2p(corr)

# specifying  Copula's
# elliptical
  if(elliptal_copula == "t") {
    # warning
    if(is.null(df_ellip))stop('Please supply a valid degrees of freedom parameter when using elliptal_copula = "t". ')
    
    Ecop <- ellipCopula(family = elliptal_copula, dispstr = "un", df = df_ellip,
                        param = Cor, dim = N)
    
    }else
      if(elliptal_copula == "norm"){
        
        Ecop <- ellipCopula(family = "norm", dispstr = "un", param = Cor, dim = N)
      
    }else
      stop("Please supply a valid argument for elliptal_copula")

# Left-cop (Archemedian copula)
  Acop <- archmCopula(family = "clayton", param = left_cop_param, dim = N)
  
  
#generating random draws from copula's, uniformly distributed data
  if(left_cop_weight<0|left_cop_weight>1)stop("Please provide a valid left_cop_weight between 0 and 1")   
  
  if(left_cop_weight==0) {
    data <- rCopula(T, Ecop)
 } else
   if(left_cop_weight==1) {
     data <- rCopula(T, Acop)
 } else
   data <- left_cop_weight*rCopula(T, Acop) + (1-left_cop_weight)*rCopula(T, Ecop)
   
  #naming and converting data to tibble
  colnames(data) <- glue::glue("Asset_{1:ncol(data)}")
  data <- as_tibble(data)
    
#Converting Uniform marginal distributions to norm or sgt. 
 if(marginal_dist=="norm"){
            
             if(is.null(sd_md))stop('Please supply a valid sd_md parameter when using marginal_dist="norm".')
             data <- data %>% map_df(~qnorm(.x, sd_md))
            
            }else
              if(marginal_dist=="sgt"){
                
                if(is.null(marginal_dist_model))stop('Please supply a valid marginal_dist_model when using marginal_dist="sgt".')
                   mu <- marginal_dist_model$mu
                   mu <- marginal_dist_model$sigma
                   lambda <- marginal_dist_model$lambda
                   p <- marginal_dist_model$p
                   q <- marginal_dist_model$q
                  
                   data <- data %>% map_df(~qsgt(.x, mu =  mu, sigma = sigma, lambda = lambda, p = p, q = q))
                  
              }
              
              return(data)
            
}
```

### Testing sim\_inno

``` r
# Sourcing function, loading data and setting seed
source("code/sim_inno.R")
source("code/gen_corr.R")
load("data/sgt_low_vol.Rda")
set.seed(872154)

Corr <- gen_corr(N = 50, Clusters = "overlapping", Num_Layers = 3, Num_Clusters = c(10,5,2))
# ----------------------------------------
# Simulating data with marginal_dist="sgt"
# ----------------------------------------
sgt_pars <- as.list(sgt_low_vol$estimate)
data_sgt <- sim_inno(corr = Corr, 
                 elliptal_copula = "t",
                 df_ellip = 4,
                 left_cop_param = 10,
                 left_cop_weight = 0,
                 marginal_dist = "sgt",
                 marginal_dist_model = sgt_pars,
                 T = 500)
colnames(data_sgt) <- glue::glue("V{1:ncol(data_sgt)}")

# ----------------------------------------
# Simulating data with marginal_dist="norm"
# ----------------------------------------
data_norm <- sim_inno(corr = Corr, 
                 elliptal_copula = "t",
                 df_ellip = 4,
                 left_cop_param = 10,
                 left_cop_weight = 0,
                 marginal_dist = "norm",
                 sd_md = 0.02311859,
                 T = 500)
colnames(data_norm) <- glue::glue("V{1:ncol(data_norm)}")

# ------------------------
# Plotting Simulated Data
# ------------------------

# First: Note how the correlation matrix has been maintained
data_sgt %>% fitHeavyTail::fit_mvt() %>% .$cov %>% cov2cor() %>% ggcorrplot::ggcorrplot()
```

<img src="README_files/figure-gfm/sim_inno test-1.png" width="80%" height="80%" style="display: block; margin: auto auto auto 0;" />

``` r
# Plotting SGT data
data_sgt %>% gather(key, value, -V1) %>% 
    mutate(Cluster = case_when( key %in% c("V2","V3","V4","V5") ~ "First Cluster",
                                key %in% c("V6","V7","V8","V9","V10") ~ "Second Cluster",
                                key %in% c("V11","V12","V13","V14","V15") ~ "Outside Cluster")) %>% 
    na.omit() %>% 
    ggplot(aes(x = V1, y = value)) +
    geom_point(alpha=0.4) +
    geom_density_2d_filled(alpha = 0.7, bins = 10) +
    facet_wrap(~Cluster, nrow = 3, 
               scales = "free_y") +
    theme_bw() +
    labs(title = "2D Density plot of SGT Distributed Data")
```

<img src="README_files/figure-gfm/sim_inno test-2.png" width="80%" height="80%" style="display: block; margin: auto auto auto 0;" />

``` r
# Plotting Normally distributed data
data_norm %>% gather(key, value, -V1) %>% 
    mutate(Cluster = case_when( key %in% c("V2","V3","V4","V5") ~ "First Cluster",
                                key %in% c("V6","V7","V8","V9","V10") ~ "Second Cluster",
                                key %in% c("V11","V12","V13","V14","V15") ~ "Outside Cluster")) %>% 
    na.omit() %>% 
    ggplot(aes(x = V1, y = value)) +
    geom_point(alpha=0.4) +
    geom_density_2d_filled(alpha = 0.7, bins = 10) +
    facet_wrap(~Cluster, nrow = 3, 
               scales = "free_y") +
    theme_bw() +
    labs(title = "2D Density plot of Normally Distributed Data")
```

<img src="README_files/figure-gfm/sim_inno test-3.png" width="80%" height="80%" style="display: block; margin: auto auto auto 0;" />

# Step 3: Introducing Volitility Persistence

The simulated innovations do not yet demonstrate the mean and/or
volatility persistence observed in real asset return series, hence why I
refer to them as innovations.

``` r
p1 <- ggAcf(data_sgt$V1) + theme_bw() + labs(title = "ACF of Innovations")
p2 <- ggAcf(data_sgt$V1^2) + theme_bw() + labs(title = "ACF of Squared Innovations")
p <- p1 / p2 
p + plot_annotation(title = "No Significant Persistence in Mean or Volatility")
```

<img src="README_files/figure-gfm/archlm-1.png" width="80%" height="80%" />

In this step I introduce autocorrelation and volatility using an AR(p,q)
+ APGARCH(q,p) model.

-   “The leptokurtosis, clustering volatility and leverage effects
    characteristics of financial time series justifies the GARCH
    modelling approach. The non-linear characteristic of the time series
    is used to check the Brownian motion and investigate into the
    temporal evolutionary patterns. The nonlinear methods of forecasting
    and signal analysis are gaining popularity in stock market because
    of their robustness in feature extraction and classiﬁcation.”
    source:
    <https://towardsdatascience.com/garch-processes-monte-carlo-simulations-for-analytical-forecast-27edf77b2787>

## sim\_garch

This function introduces mean and variance persistence by plugging in
the numbers generated by the sim\_inno function as innovations in the
GARCH process. Note that most of code was borrowed from fGarch’s
garchspec and garchsim functions.

``` r
sim_garch <- function(model= list(), innovations, simple = TRUE){
    
    #default parameters for garch model  
    default <- list(omega = 1e-06, 
                    alpha = 0.1,
                    gamma = NULL, 
                    beta = 0.8, 
                    mu = 0,   #changed form NULL to 0
                    ar = NULL, 
                    ma = NULL, 
                    delta = 2)
    
    default[names(model)] <- model
    model <- default   
    
    #obtaining parameters and lag orders from model object
    mu <- model$mu 
    ar <- model$ar
    ma <- model$ma
    omega <-  model$omega
    alpha <- model$alpha
    gamma <- model$gamma
    beta <- model$beta
    delta <- model$delta
    deltainv <- 1/delta
    order.ar <- length(ar)
    order.ma <- length(ma)
    order.alpha <- length(alpha)
    order.beta <- length(beta)
    max.order <- max(order.ar, order.ma, order.alpha, order.beta)
    n <- length(innovations) - 5
    
    if(max.order>5)stop("Please supply a volitility model with max order less than or equal to 5")
    
    #Generating innovations
    z_length <- n + max.order
    z <- c(innovations)[1:z_length]  #must change this later
    
    h <- c(rep(model$omega/(1 - sum(model$alpha) - sum(model$beta)), 
               times = max.order), rep(NA, n))    #sd's
    
    y <- c(rep(model$mu/(1 - sum(model$ar)), times = max.order), rep(NA, n))  #garch simulations
    m <- max.order
    
    #simulating GARCH  
    eps <- h^deltainv * z  #this part often breaks depending on GARCH parameters chosen
    
    for (i in (m + 1):(n + m)) {
        h[i] = omega + sum(alpha * (abs(eps[i - (1:order.alpha)]) - 
                                        gamma * (eps[i - (1:order.alpha)]))^delta) + 
            sum(beta * h[i - (1:order.beta)])  
        
        eps[i] = h[i]^deltainv * z[i]
        y[i] = mu + sum(ar * y[i - (1:order.ar)]) + sum(ma * 
                                                            eps[i - (1:order.ma)]) + eps[i]
    }
    
    if(simple == TRUE) {
        data <- y[(m + 1):(n + m)] #removes burn in data
    } else {
        data <- tibble(z = z[(m + 1):(n + m)], 
                       sigma = h[(m + 1):(n + m)]^deltainv, 
                       y = y[(m + 1):(n + m)])
    }
    
}
```

## Demonstrating sim\_garch

Note that: - volatility clusters and significant autocorrelation appear
in the series after it is processed through the sim\_garch function. -
sim\_garch is a deterministic function.

``` r
source("code/sim_inno.R")
load("data/garch.Rda")
set.seed(32156454)

inno <- sgt::rsgt(n = 10001, lambda = -0.0143, p = 1.6650, q = 1.9095)
#Parameters from Statistics and Data Analysis for Financial Engineering pg.421-423
model <- list(mu = 0.000002,        
              omega = 0.000005, #key unconditional volatility parameter
              alpha = 0.098839, 
              beta = 0.899506, 
              ar = 0.063666,
              ma = NULL,
              gamma = 0.12194,
              delta = 1.85)
# Introducing Mean and Var Persistence
return <- sim_garch(model, inno, simple = TRUE)

p_load(patchwork)
p1 <- inno %>% as_tibble() %>% ggplot(aes(x=1:length(inno), y=value)) +
  geom_line() + theme_bw() + labs(subtitle = "Random Draws From SGT Distribution", x = "", y = "Innovations")

p2 <- return %>% as_tibble() %>% ggplot(aes(x=1:length(return), y=value)) +
  geom_line() + theme_bw() + labs(subtitle = "Same Random Draws After sim_garch", x = "", y = "Returns")
  
p1/p2 + plot_annotation(title = "SGT Innovations vs APGARCH Returns")
```

<img src="README_files/figure-gfm/unnamed-chunk-4-1.png" width="80%" height="80%" />

``` r
p1 <- ggAcf(inno^2)+ theme_bw()+ labs(title = "ACF of Squared Innovations")
p2 <- ggAcf(return^2) + theme_bw() + labs(title = "ACF of Squared Returns")
p1/p2
```

<img src="README_files/figure-gfm/unnamed-chunk-4-2.png" width="80%" height="80%" />

# Steps 1, 2 and 3: Simulating an Asset Market

The code below uses the sim\_inno and sim\_garch functions to simulate
500 days of returns for a market of 20 assets.

``` r
source("code/gen_corr.R")
source("code/sim_inno.R")

set.seed(123)

# toy corr matrix
# corr <- gen_corr(N = 20, Clusters = "none", Num_Layers = 3, Num_Clusters = c(2,4,5))   

# Empirical Correlation matrix
load("data/labeled_training_data.Rda")
set.seed(1234)
corr <- labeled_training_data$stressed_market[[1]]
dim <- sample(1:nrow(corr), size = 20)
corr <- corr[dim,dim]

# marginal distribution parameters (for SGT)
sgt_pars <- list(sigma = 1, lambda = -0.04140381, p = 1.880649, q = 1.621578)

inno <- sim_inno(corr = corr, 
         elliptal_copula = "t",
         df_ellip = 4,
         left_cop_param = 2,
         left_cop_weight = 0.01,
         marginal_dist = "sgt",
         marginal_dist_model = sgt_pars,
         T = 500)

#getting empirical parameters
load("data/coef.Rda")

model <- coef %>% select_if(is.numeric) %>% map_df(mean)

#Parameters from Statistics and Data Analysis for Financial Engineering pg.421-423
model <- list(mu = 0.000002,        
              omega = 0.000005, #key unconditional volatility parameter
              alpha = 0.098839, 
              beta = 0.899506, 
              ar = 0.063666,
              ma = NULL,
              gamma = 0.12194,
              delta = 1.85)


#Making Simdat tidy
simdat <- inno %>% map_dfc(~sim_garch(model, .x))

tidy_simdat <- simdat %>%  mutate(date = 1:nrow(simdat)) %>% 
  gather(key = Asset, value = Return, -date)

#Correlation matrix of GARCH sim data
#simdat %>% spread(key = Asset, value = Return) %>% select(-date) %>% cor() %>% corrplot::corrplot()

#generating and plotting Cum returns
tidy_simdat %>%  
  arrange(date) %>% 
  group_by(Asset) %>%
  mutate(Cum_Return = cumprod(1 + Return)*100) %>%
  ggplot() + geom_line(aes(x=date, y=Cum_Return, color = Asset)) +
  facet_wrap(~Asset) +
  labs(title = "Cumulative Returns",
       subtitle = "500 Trading Days Across 20 Assets") +
  theme_bw() +
  theme(legend.position = "none") 
```

<img src="README_files/figure-gfm/sim market-1.png" width="80%" height="80%" />

``` r
tidy_simdat %>% ggplot() +
  geom_line(aes(x=date,y=Return, color = Asset)) +
  facet_wrap(~Asset) +
  labs(title = "Returns",
       subtitle = "500 Trading Days Across 20 Assets") +
  theme_bw() +
  theme(legend.position = "none")
```

<img src="README_files/figure-gfm/sim market-2.png" width="80%" height="80%" />

## Prototype: sim\_asset\_market Function

``` r
source("code/gen_corr.R")
source("code/sim_inno.R")
sim_asset_market <- function(corr, T = 500, model = list()){    #note that length out is 499 not 500,  

        sgt_pars <- list(mu = 0, sigma = 1, 
                         lambda = -0.04140381, p = 1.880649, q = 1.621578)
        inno <- sim_inno(corr = corr, 
                 elliptal_copula = "t",
                 df_ellip = 4,
                 left_cop_param = 4,
                 left_cop_weight = 0,
                 marginal_dist = "sgt",
                 marginal_dist_model = sgt_pars,
                 T = T)
        
        #Applying sim_garch to each column in simdat
        simdat <- inno %>% map_dfc(~sim_garch(model, .x))
        
        #adding a date column
       start_date <- Sys.Date()
       all_days <- seq(start_date, start_date %m+% lubridate::days( ceiling( T*(1+(3/7))) ), by = 1)
       weekdays <- all_days[!weekdays(all_days) %in% c('Saturday','Sunday')][1:T]
        
        #Creating final df
        simdat <- simdat %>% 
            mutate(date = weekdays, .before = `Asset_1`) %>% 
            gather(key=Asset, value = Return, -date)
        
    return(simdat)
}
```

# Step 4: Simulating an Ensemble of Asset Markets

## Now lets use the sim\_asset\_market above to the first MC simulation.

``` r
#Emperical Corr matrix
load("data/labeled_training_data.Rda")

set.seed(2121354)
corr <- labeled_training_data$stressed_market[[1]]
dim <- sample(1:nrow(corr), size = 20)
corr <- corr[dim,dim]
ggcorrplot::ggcorrplot(corr, hc.order = TRUE)
```

<img src="README_files/figure-gfm/monte carlo-1.png" width="80%" height="80%" />

``` r
#Parameters from Statistics and Data Analysis for Financial Engineering pg.421-423
model <- list(mu = 0,        
              omega = 0.0000025, #key unconditional volatility parameter
              alpha = 0.098839, 
              beta = 0.899506, 
              ar = 0.063666,
              ma = NULL,
              gamma = 0.12194,
              delta = 1.85)

# Simulating N markets
N <- 50
# How can I speed this up?
mc_data <- 1:N %>% map(~sim_asset_market(corr, T = 500, model)) %>% 
  reduce(left_join, by = c("date","Asset")) 

# Setting appropriate column names
colnames(mc_data) <- c("date", "Asset", glue::glue("Universe_{1:(ncol(mc_data)-2)}"))

# making data tidy: This is how I want my final output to look. 
mc_data <- mc_data %>% gather(Universe, Return, c(-date, -Asset))

# Plotting reasults
mc_data %>% 
  group_by(Asset, Universe) %>% 
  arrange(date) %>% 
  mutate(cum_ret = cumprod(1 + Return)*100) %>% 
  ggplot() +
  geom_line(aes(x = date, y = cum_ret, color = Universe), size = 1, alpha = 0.5) + 
  facet_wrap(~Asset, scales = "free_y") + 
  labs(title = "Ensemble of Cumulative Returns",
       subtitle = "50 Realizations for a Market of 20 Assets") +
  theme_bw()+
  theme(legend.position = "none") 
```

<img src="README_files/figure-gfm/monte carlo-2.png" width="80%" height="80%" />

## mc\_market

``` r
mc_market <- function(){
  
}
```
