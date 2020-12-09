sim_inno <- function(corr,
                     T = 252,
                     mv_dist = "t",
                     df_ellip = 3,
                     left_cop_param = 5,
                     left_cop_weight = 0,
                     marginal_dist = NULL,
                     marginal_dist_model = NULL) {

    N <- nrow(corr)
    T <- T + 5   # extra room for GARCH to use later
    Cor <- P2p(corr)

    # Specifying  Copulas
    # elliptical
    if (mv_dist == "t") {
        # warning
        if (is.null(df_ellip)) stop('Please supply a valid degrees of freedom parameter when using mv_dist = "t".')

        Ecop <- ellipCopula(family = mv_dist, dispstr = "un", df = df_ellip,
                            param = Cor, dim = N)

    } else
        if (mv_dist == "norm") {

            Ecop <- ellipCopula(family = "norm", dispstr = "un", param = Cor, dim = N)

        } else
            if(!(mv_dist %in% c("norm", "t"))) stop("Please supply a valid argument for mv_dist")

    # Left-cop (Archemedian copula)
    if (left_cop_weight != 0) {
        Acop <- archmCopula(family = "clayton", param = left_cop_param, dim = N)
    }

    #generating random (uniformly distributed) draws from hybrid copula's
    if (left_cop_weight<0|left_cop_weight>1) stop("Please provide a valid left_cop_weight between 0 and 1")

    if (left_cop_weight==0) {
        data <- rCopula(T, Ecop)
    } else
        if(left_cop_weight==1) {
            data <- rCopula(T, Acop)
        } else
            data <- left_cop_weight*rCopula(T, Acop) + (1-left_cop_weight)*rCopula(T, Ecop)

    #naming and converting data to tibble
    colnames(data) <- glue::glue("Asset_{1:ncol(data)}")
    data <- as_tibble(data)

    if (is.null(marginal_dist)) return(data)

    #Converting Uniform marginal distributions to norm or sgt.
    if (!(marginal_dist %in% c("norm", "t", "sgt"))) stop ("Please supply a valid marginal_dist argument")
    if (marginal_dist=="norm") {
        if (is.null(marginal_dist_model)) {
            marginal_dist_model <- list(mu=0, sigma = 1)
        }
        mu <- marginal_dist_model$mu
        sigma <- marginal_dist_model$sigma
        data <- data %>% map_df(~qnorm(.x, mean = mu, sd = sigma))
        return(data)
    } else
        if (marginal_dist == "t"){
            mu <- marginal_dist_model$mu
            df <- marginal_dist_model$df
            data <- data %>%
                map_df(~qt(.x, ncp = mu, df = df))
            return(data)
        } else
            if (marginal_dist == "sgt") {
                if (is.null(marginal_dist_model))
                    stop ('Please supply a valid marginal_dist_model when using marginal_dist="sgt".')
                mu <- marginal_dist_model$mu
                sigma <- marginal_dist_model$sigma
                lambda <- marginal_dist_model$lambda
                p <- marginal_dist_model$p
                q <- marginal_dist_model$q
                data <- data %>%
                    map_df(~qsgt(.x, mu =  mu, sigma = sigma, lambda = lambda, p = p, q = q, mean.cent = TRUE))
                return(data)
                }
}
