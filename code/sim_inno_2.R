sim_inno_2 <- function(corr,
                     k = 252,
                     mv_dist = "t",
                     mv_df = 4,
                     left_cop_weight = 0,
                     left_cop_param = 5,
                     marginal_dist = "norm",
                     marginal_dist_model = NULL,
                     ts_model = NULL) {

    N <- nrow(corr)
    k <- k + 5   # extra room for sim_garch to as lags at later stage.
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

        data <- data %>% dplyr::filter(date >= first(date) %m+% days(5))
        return(data)

    } else {
        data <- data %>% group_by(Asset) %>% arrange(date) %>%
            mutate(Return = sim_garch(Return, model = ts_model)) %>% na.omit()
        return(data)
        }
}

set.seed(123)
sim_inno_2(diag(10), k = 300, mv_dist = "norm", mv_df = 2, left_cop_weight = 0.5,
           left_cop_param = 2, marginal_dist = "sgt",
           marginal_dist_model = list(mu = 0, sd = 1, lambda = 0.5, p = Inf, q = 2),
           ts_model = list(omega = 1))
set.seed(123)
sim_inno_2(diag(10), k = 300, mv_dist = "norm", mv_df = 2, left_cop_weight = 0.5,
           left_cop_param = 2, marginal_dist = "sgt",
           marginal_dist_model = list(mu = 0, sd = 1, lambda = 0.5, p = Inf, q = 2),
           ts_model = NULL)




#Edit sim garch to give NA's for burn in period: Done!!!



