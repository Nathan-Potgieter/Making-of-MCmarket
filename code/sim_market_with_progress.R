#' @title sim_market_with_progress
#' @description This is an alternative version of the sim_market function that
#' includes functionality to include a progress bar when used in purrr::map functions
#' @note  See ??sim_market for details on use.
#' @return a tidy tibble containing a date, Asset and Return column.
#'
#' @importFrom copula ellipCopula archmCopula rcopula
#' @importFrom glue glue
#' @importFrom sgt qsgt
#' @importFrom lubridate '%m+%' days
#' @importFrom tidyr gather
#' @import dplyr
#' @import purrr
#'
#' @examples
#'
#' \dontrun{
#'
#' library(MCmarket)
#' library(tidyverse)
#'
#' ### creating a correlation matrix to use as input in sim_asset_market
#' corr <- gen_corr(N = 20, Clusters = "none")
#'
#'
#' N <- 100
#' pb <- dplyr::progress_estimated(N)   # this must be named pb
#' market <-
#'       map_dfr(1:N,
#'               ~sim_market_with_progress(corr),
#'               .id = "Universe")
#'
#' }
#' @export
#'
sim_market_with_progress <- function(corr,
                                     k = 252,
                                     mv_dist = "t",
                                     mv_df = 3,
                                     left_cop_weight = 0,
                                     left_cop_param = 4,
                                     marginal_dist = "norm",
                                     marginal_dist_model = NULL,
                                     ts_model = list()
) {

    # Tick to progress bar
    pb$tick()$print()

    # sim_market
    N <- nrow(corr)
    k <- k + 1   # extra room for sim_garch to as a lag.
    Cor <- P2p(corr)

    # Specifying  Copulas
    # elliptical
    if(!(mv_dist %in% c("norm", "t"))) stop("Please supply a valid argument for mv_dist")
    else
        if (mv_dist == "t") {
            if (is.null(mv_df)) stop('Please supply a valid degrees of freedom parameter when using mv_dist = "t".')
            Ecop <- ellipCopula(family = "t",
                                dispstr = "un",
                                df = mv_df,
                                param = Cor,
                                dim = N)
        } else
            if (mv_dist == "norm") {
                Ecop <- ellipCopula(family = "normal",
                                    dispstr = "un",
                                    param = Cor,
                                    dim = N)
            }

    # Left-cop (Archemedian copula)
    if (left_cop_weight < 0|left_cop_weight > 1) stop("Please provide a valid left_cop_weight between 0 and 1")
    if (left_cop_weight != 0) {
        Acop <- archmCopula(family = "clayton",
                            param = left_cop_param,
                            dim = N)
    }

    # Generating random (uniformly distributed) draws from hybrid copula's
    if (left_cop_weight == 0) {
        data <- rCopula(k, Ecop)
    } else
        if(left_cop_weight == 1) {
            data <- rCopula(k, Acop)
        } else {
            data <- (left_cop_weight*rCopula(k, Acop) + (1-left_cop_weight)*rCopula(k, Ecop))
        }


    # Creating a date vector
    start_date <- Sys.Date()
    dates <- rmsfuns::dateconverter(StartDate = start_date,
                                    EndDate = start_date %m+% lubridate::days(k-1),
                                    Transform = "alldays")

    # Making Tidy & adding date column
    data <- as_tibble(data) %>%
        purrr::set_names(glue::glue("Asset_{1:ncol(data)}")) %>%
        mutate(date = dates) %>%
        gather(Asset, Value, -date)


    if (!(marginal_dist %in% c("norm", "t", "sgt", "unif"))) stop ("Please supply a valid marginal_dist argument")

    if (marginal_dist == "unif") return(data)

    # Warnings
    if (marginal_dist == "norm" & is.null(marginal_dist_model)) marginal_dist_model <- list(mu=0, sd = 1)
    if (marginal_dist == "t" & is.null(marginal_dist_model))  marginal_dist_model <- list(mu=0, df = 5)
    if (marginal_dist == "sgt" & is.null(marginal_dist_model)) stop ('Please supply a valid marginal_dist_model when using marginal_dist="sgt".')

    #Converting Uniform marginal distributions to norm, t or sgt.
    args <- tibble(Asset = glue::glue("Asset_{1:N}")) %>%
        mutate(mean = marginal_dist_model$mu,
               sd = marginal_dist_model$sd,
               ncp = marginal_dist_model$ncp,
               df = marginal_dist_model$df,
               lambda = marginal_dist_model$lambda,
               p = marginal_dist_model$p,
               q = marginal_dist_model$q)

    if (marginal_dist == "norm") {

        if(is.null(marginal_dist_model$mu)) stop('Please supply a valid mu parameter when using marginal_dist = "norm".')
        if(is.null(marginal_dist_model$sd)) stop('Please supply a valid sd parameter when using marginal_dist = "norm".')

        data <- data %>% left_join(., args, by = "Asset") %>%
            group_by(Asset) %>%  arrange(date) %>%
            mutate(Return =  qnorm(Value, mean, sd)) %>%
            select(date, Asset, Return)

    } else
        if (marginal_dist == "t") {

            if(is.null(marginal_dist_model$ncp)) stop('Please supply a valid ncp parameter when using marginal_dist = "t".')
            if(is.null(marginal_dist_model$df)) stop('Please supply a valid df parameter when using marginal_dist = "t".')

            data <- data %>% left_join(., args, by = "Asset") %>%
                group_by(Asset) %>%  arrange(date) %>%
                mutate(Return = qt(Value, df =  df, ncp =  ncp)) %>%
                select(date, Asset, Return)

        } else
            if (marginal_dist == "sgt") {

                if (is.null(marginal_dist_model$mu)) marginal_dist_model$mu <- 0
                if (is.null(marginal_dist_model$sd)) marginal_dist_model$sd <- 1
                if (is.null(marginal_dist_model$lambda)|
                    is.null(marginal_dist_model$p)|
                    is.null(marginal_dist_model$q)) stop('Please supply valid arguments for lambda, p and q when using marginal_dist = "sgt".')

                data <- data %>% left_join(., args, by = "Asset") %>%
                    group_by(Asset) %>% arrange(date) %>%
                    mutate(Return = qsgt(Value, mean, sd, lambda, p, q)) %>%
                    select(date, Asset, Return)

            }

    if (is.null(ts_model)) {

        data <- data %>% dplyr::filter(date > first(date))
        return(data)

    } else {

        # Warnings
        if (!is.null(ts_model$omega) & length(ts_model$omega) != 1 & length(ts_model$omega) != N) stop("Please supply a valid vector length for omega. Must be of length 1 or ncol(corr).")
        if (!is.null(ts_model$alpha) & length(ts_model$alpha) != 1 & length(ts_model$alpha) != N) stop("Please supply a valid vector length for alpha. Must be of length 1 or ncol(corr).")
        if (!is.null(ts_model$gamma) & length(ts_model$gamma) != 1 & length(ts_model$gamma) != N) stop("Please supply a valid vector length for gamma. Must be of length 1 or ncol(corr).")
        if (!is.null(ts_model$beta) & length(ts_model$beta) != 1 & length(ts_model$beta) != N) stop("Please supply a valid vector length for beta. Must be of length 1 or ncol(corr).")
        if (!is.null(ts_model$mu) & length(ts_model$mu) != 1 & length(ts_model$mu) != N) stop("Please supply a valid vector length for mu. Must be of length 1 or ncol(corr).")
        if (!is.null(ts_model$ar) & length(ts_model$ar) != 1 & length(ts_model$ar) != N) stop("Please supply a valid vector length for ar. Must be of length 1 or ncol(corr).")
        if (!is.null(ts_model$ma) & length(ts_model$ma) != 1 & length(ts_model$ma) != N) stop("Please supply a valid vector length for ma. Must be of length 1 or ncol(corr).")
        if (!is.null(ts_model$delta) & length(ts_model$delta) != 1 & length(ts_model$delta) != N) stop("Please supply a valid vector length for delta. Must be of length 1 or ncol(corr).")

        # Inducing mean and/or variance persistence

        # Tibble with with garh parameters and defaults
        ts_args <- tibble(Asset = glue::glue("Asset_{1:N}")) %>%
            mutate(omega = if (is.null(ts_model$omega)) {5e-04} else {ts_model$omega},
                   alpha = if (is.null(ts_model$alpha)) {0} else {ts_model$alpha},
                   gamma = if (is.null(ts_model$gamma)) {0}  else {ts_model$gamma},
                   beta = if (is.null(ts_model$beta)) {0} else {ts_model$beta},
                   mu = if (is.null(ts_model$mu)) {0} else {ts_model$mu},   #changed form NULL to 0
                   ar = if (is.null(ts_model$ar)) {0} else {ts_model$ar},
                   ma = if (is.null(ts_model$ma)) {0} else {ts_model$ma},
                   delta = if (is.null(ts_model$delta)) {2} else {ts_model$delta})

        # Inducing garch properties
        data <- data %>% left_join(., ts_args, by = "Asset") %>%
            arrange(date) %>% group_by(Asset) %>%
            mutate(Return = sim_garch(innovations = Return,
                                      omega = omega,
                                      alpha = alpha,
                                      gamma = gamma,
                                      beta = beta,
                                      mu = mu,
                                      ar = ar,
                                      ma = ma,
                                      delta = delta)) %>% na.omit() %>%
            select(date, Asset, Return)
        return(data)
    }
}
