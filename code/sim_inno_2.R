sim_inno_2 <- function(corr,
                     k = 252,
                     mv_dist = "t",
                     mv_df = 4,
                     left_cop_weight = 0,
                     left_cop_param = 5,
                     marginal_dist = "norm",
                     marginal_dist_model = NULL) {

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
        data <- rCopula(k, Ecop) %>% as_tibble() %>%
            purrr::set_names(glue::glue("Asset_{1:ncol(data)}")) %>%
            gather(Asset, Return)
    } else
        if(left_cop_weight == 1) {
            data <- rCopula(k, Acop) %>% as_tibble() %>%
                purrr::set_names(glue::glue("Asset_{1:ncol(data)}")) %>%
                gather(Asset, Return)
        } else
            data <- (left_cop_weight*rCopula(k, Acop) + (1-left_cop_weight)*rCopula(k, Ecop)) %>%
        as_tibble() %>% purrr::set_names(glue::glue("Asset_{1:ncol(data)}")) %>%
        gather(Asset, Return)


    if (!(marginal_dist %in% c("norm", "t", "sgt", "unif"))) stop ("Please supply a valid marginal_dist argument")

    if (marginal_dist == "unif") return(data)

    #Converting Uniform marginal distributions to norm, t or sgt.

    args <- data.frame(Type = glue::glue("Asset_{1:ncol(data)}"),
                       mean = marginal_dist_model$mu,
                       sd = marginal_dist_model$sd)





    if (marginal_dist == "norm") {
        if (is.null(marginal_dist_model)) {
            marginal_dist_model <- list(mu=0, sd = 1)
        }

        data <- data %>% map_df(~qnorm(.x,
                                       mean = marginal_dist_model$mu,
                                       sd = marginal_dist_model$sd))
        return(data)

    } else
        if (marginal_dist == "t") {
            if (is.null(marginal_dist_model)) {
                marginal_dist_model <- list(mu=0, df = 5)
            }

            data <- data %>%
                map_df(~qt(.x,
                           ncp = marginal_dist_model$mu,
                           df = mv_df))
            return(data)

        } else
            if (marginal_dist == "sgt") {
                if (is.null(marginal_dist_model))
                    stop ('Please supply a valid marginal_dist_model when using marginal_dist="sgt".')
                else
                    if (is.null(marginal_dist_model$mu)) {marginal_dist_model$mu <- 0}
                if (is.null(marginal_dist_model$sd)) {marginal_dist_model$sd <- 1}  # Do I need an else?? Seems like including it will not work
                if (is.null(marginal_dist_model$lambda)|
                    is.null(marginal_dist_model$p)|
                    is.null(marginal_dist_model$q)) stop('Please supply valid arguments for lambda, p and q when using marginal_dist = "sgt".')

                data <- data %>%
                    map_df(~qsgt(.x,
                                 mu =  marginal_dist_model$mu,
                                 sigma = marginal_dist_model$sd,
                                 lambda = marginal_dist_model$lambda,
                                 p = marginal_dist_model$p,
                                 q = marginal_dist_model$q,
                                 mean.cent = TRUE))
                return(data)
            }
}
