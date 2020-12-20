#' @title sim_inno
#' @description This function simulates a multivariate data set with a given correlation structure. The multivariate distribution can be normal or student t distributed, or a weighted average between one and the Clayton copula's multivariate distribution. The univariate distributions can set to normal, student t, skewed generalized t or uniform distributions.
#' Note the function name sim_inno, simulate innovations, stems from the i.i.d nature of the resulting series. Use the sim_asset_market function if you require the data to exhibit time-series persistence.
#' @param corr a correlation matrix that the simulated date will adhere to. Note that the number of variables simulated is equal to the number of columns in the correlation matrix.
#' @param k a positive integer indicating the number of time periods to simulate. Note that the number of periods generated is actually equal to k + 5 as these extra observations are needed when applying time series properties to the data.
#' @param mv_dist a string specifying the multivariate distribution. Can be one of c("norm", "t") referring to the multivariate normal and t distributions respectively. Default is 3.
#' @param mv_df degrees of freedom of the multivariate distribution, required when mv_dist = "t".
#' @param left_cop_weight a positive value between zero and one indicating the weight applied to the Clayton copula when creating the multivariate distribution. Note that a value between zero and one essentially generates a hybrid distribution between the chosen mv_dist and the Clayton copula. Therefore, the greater the left_cop_weight the less the data will reflect the correlation structure. Default is set to 0.
#' @param left_cop_param a positive value indicating the parameter of the Clayton copula. Default is 4.
#' @param marginal_dist a string variable specifying the univariate distribution of each variable. Can be one of c("norm", "t", "sgt") referring to the normal, student-t and skewed-generalized-t distributions respectively. Default is "norm".
#' @param  marginal_dist_model list containing the relevant parameters for the chosen marginal_dist. marginal_dist = "norm" accepts a mean and standard deviation with defaults list(mu = 0, sigma = 1) respectively. marginal_dist = "t" accepts the non-centrality and degrees of freedom arguments, default values are list(mu = 0, df = 5). marginal_dist = "sgt" accepts the mean, sd, lambda, p and q parameters list(mu = 0, sigma = 1, lambda, p, q). Note lambda, p and q have no defaults and must therefore be set by the user.
#' @return A tibble with columns containing the values for each variable.
#'
#' @importFrom xts as.xts timeBased
#'
#' @import copula
#' @import glue
#' @import dplyr
#' @import purrr
#' @import sgt
#'
#' @importFrom rlang :=
#' @importFrom rlang := enquo quo_get_expr
#' @examples
#'
#' \dontrun{
#' library(copula)
#' library(glue)
#' library(dplyr)
#' library(purrr)
#' library(sgt)
#'
#' ### generating a correlation matrix to be used
#' corr <- gen_corr()
#'
#' ### simulating 252 observations multivariate t distributed date with sgt univariate distributions.
#' data <- sim_inno(corr,
#'                  k = 252,
#'                  mv_dist = "t",
#'                  mv_df = "5"
#'                  marginal_dist = "sgt",
#'                  marginal_dist_model = list(mu = 0, sigma = 1,
#'                                             lambda = -0.1, p = 2 ,
#'                                             q = Inf))
#' print(data)
#'
#' }
#' @export

sim_inno <- function(corr,
                     k = 252,
                     mv_dist = "t",
                     mv_df = 3,
                     left_cop_weight = 0,
                     left_cop_param = 5,
                     marginal_dist = "norm",
                     marginal_dist_model = NULL) {

    N <- nrow(corr)
    k <- k + 5   # extra room for sim_garch to use later
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
            Ecop <- ellipCopula(family = "normal", dispstr = "un", param = Cor, dim = N)
        }

    # Left-cop (Archemedian copula)
    if (left_cop_weight != 0) {
        Acop <- archmCopula(family = "clayton", param = left_cop_param, dim = N)
    }

    #generating random (uniformly distributed) draws from hybrid copula's
    if (left_cop_weight < 0|left_cop_weight > 1) stop("Please provide a valid left_cop_weight between 0 and 1")

    if (left_cop_weight == 0) {
        data <- rCopula(k, Ecop)
    } else
        if(left_cop_weight == 1) {
            data <- rCopula(k, Acop)
        } else
            data <- left_cop_weight*rCopula(k, Acop) + (1-left_cop_weight)*rCopula(k, Ecop)

    #naming and converting data to tibble
    colnames(data) <- glue::glue("Asset_{1:ncol(data)}")
    data <- as_tibble(data)

    if (!(marginal_dist %in% c("norm", "t", "sgt", "unif"))) stop ("Please supply a valid marginal_dist argument")

    if (marginal_dist == "unif") return(data)  # Does function stop here; Do I put an else here?

    #Converting Uniform marginal distributions to norm, t or sgt.

    if (marginal_dist == "norm") {
        if (is.null(marginal_dist_model)) {
            marginal_dist_model <- list(mu=0, sigma = 1)
        }

        data <- data %>% map_df(~qnorm(.x,
                                       mean = marginal_dist_model$mu,
                                       sd = marginal_dist_model$sigma))
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
                    if (is.null(marginal_dist_model$sigma)) {marginal_dist_model$sigma <- 1}  # Do I need an else?? Seems like including it will not work
                    if (is.null(marginal_dist_model$lambda)|
                       is.null(marginal_dist_model$p)|
                       is.null(marginal_dist_model$q)) stop('Please supply valid arguments for lambda, p and q when using marginal_dist = "sgt".')

                        data <- data %>%
                            map_df(~qsgt(.x,
                                         mu =  marginal_dist_model$mu,
                                         sigma = marginal_dist_model$sigma,
                                         lambda = marginal_dist_model$lambda,
                                         p = marginal_dist_model$p,
                                         q = marginal_dist_model$q,
                                         mean.cent = TRUE))
                        return(data)
                }
}
