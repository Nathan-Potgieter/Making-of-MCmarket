
gen_inno <- function(corr,
                     elliptal_copula = c("norm", "t"),
                     df_ellip = NULL,
                     left_cop_param = 5,
                     left_cop_weight = 0.5,
                     T = 251,
                     marginal_dist = NULL,
                     sd_md = NULL,
                     nu_md = NULL,
                     lambda_md = 0,
                     p_md = 2,
                     q_md = Inf){

    N <- nrow(corr)
    Cor <- P2p(corr)

    #specifying  Copula's
    #elliptal
    if(elliptal_copula == "t"){

        if(is.null(df_ellip))stop('Please supply a valid degrees of freedom parameter when using elliptal_copula = "t". ')
        Ecop <- ellipCopula(family = elliptal_copula, dispstr = "un", df = df_ellip, param = Cor, dim = N)

    }else
        if(elliptal_copula == "norm"){

            Ecop <- ellipCopula(family = "norm", dispstr = "un", param = Cor, dim = N)

        }else stop("Please supply a valid argument for elliptal_copula")


    #left-cop

    Acop <- archmCopula(family = "clayton", param = left_cop_param, dim = N)





    data <- left_cop_weight*rCopula(T, Acop) + (1-left_cop_weight)*rCopula(T, Ecop)


    #Converting Uniform marginal distributions to t or norm
    if(marginal_dist=="std"){

        if(is.null(sd_md)|is.null(nu_md))stop('Please supply a valid sd and nu parameter when using marginal_dist=="std". ')
        data <- apply(data, 2, qstd, sd = sd_md, nu = nu_md)

    }else
        if(marginal_dist=="norm"){

            if(is.null(sd_md))stop('Please supply a valid sd parameter when using marginal_dist="norm".')
            data <- apply(data, 2, qnorm, sd = sd_md)

        }else
            if(marginal_dist=="sgt"){

                if(is.null(sd_md))stop('Please supply a valid sd_md parameter when using marginal_dist="sgt".')
                data <- apply(data, 2, qsgt, sigma = sd_md, lambda = lambda_md, p = p_md, q = q_md)

            }

    return(data)

}