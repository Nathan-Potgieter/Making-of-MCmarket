hycop <- function(corr,
                  elliptal_copula = c("norm", "t"),
                  df_ellip = NULL,
                  left_cop_param = 5,
                  left_cop_weight = 0.5,
                  T = 251,
                  marginal_dist = NULL,
                  df_marginal_dist = NULL){

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
    if(is.null(marginal_dist)==TRUE){
        return(data)
    }else
        if(marginal_dist=="t"){

            if(is.null(df_marginal_dist))stop('Please supply a valid degrees of freedom parameter (df_marginal_dist) when using marginal_dist=="t". ')

            data <- apply(data, 2, qt, df = df_marginal_dist)

        }else
            if(marginal_dist=="norm"){
                data <- apply(data, 2, qnorm)
            }

}