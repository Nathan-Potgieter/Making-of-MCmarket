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

            lambda <- marginal_dist_model$lambda
            p <- marginal_dist_model$p
            q <- marginal_dist_model$q

            data <- data %>% map_df(~qsgt(.x, sigma = 1, lambda = lambda, p = p, q = q))

        }

    return(data)

}