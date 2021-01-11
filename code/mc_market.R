#' @title mc_market
#' @description This function produces an ensemble of market returns with a given correlation matrix.
#' The user can choose between the multivariate t and normal
#' distributions and adjust the markets left tail dependency by weighting in the Clayton copula.
#' The univariate asset return distributions can also be set to normal, student-t or sgt
#' distributed. Finally, mean and variance persistence can be induced via the parameters of an
#' ARMA + APARCH model.
#' @note See examples under sim_market_with_progress for instructions on how to add an on screen
#'  progress bar when performing the Monte Carlo simulation, this is recommended as calculations may
#'  take a number of minuets.
#'
#' It is suggested that the marginal distributions be set to mean zero and standard deviation
#' one. Those attributes are better set in the ts_model argument.
#' @param corr a correlation matrix specifying the correlation structure of the simulated data.
#' Note that the number of variables simulated is equal to the number of columns in the correlation matrix.
#' @param N a positive integer indicating the number of markets to simulate.
#' @param k a positive integer indicating the number of time periods to simulate.
#' @param mv_dist a string specifying the multivariate distribution. Can be one of c("norm", "t")
#' which correspond to the multivariate normal and t distributions respectively.
#' @param mv_df degrees of freedom of the multivariate distribution, required when mv_dist = "t". Default is 3.
#' @param left_cop_weight a positive value between zero and one stipulating the weight applied to
#' the Clayton copula when simulating the multivariate distribution. Note that a value between zero
#' and one essentially generates a hybrid distribution between the chosen mv_dist and the Clayton
#' copula. Therefore, the greater the left_cop_weight the less the data will reflect the correlation
#' structure. Default is set to 0.
#' @param left_cop_param a positive value indicating the parameter of the Clayton copula. Default is 4.
#' @param marginal_dist a string variable specifying the asset returns univariate distribution. Can
#' be one of c("norm", "t", "sgt") referring to the normal, student-t and skewed-generalized-t
#' distributions respectively. Default is "norm".
#' @param  marginal_dist_model list containing the relevant parameters for the chosen marginal_dist.
#'
#' marginal_dist = "norm" accepts the mean (mu) and standard deviation (sd) arguments with their respective
#' defaults set to list(mu = 0, sd = 1).
#'
#' marginal_dist = "t" accepts the non-centrality parameter (ncp) and degrees of freedom (df) arguments,
#' default values are list(ncp = 0, df = 5).
#'
#' marginal_dist = "sgt" accepts the mean (mu), standard deviation (sd), lambda, p and q parameters
#' list(mu = 0, sigma = 1, lambda, p, q). Note lambda, p and q have no defaults and must therefore be set
#' by the user. For more information see the documentation on the qsgt function from the sgt pacakge.
#' @param ts_model a list containing various ARMA + APGARCH model parameters allowing one to specify
#' the time series properties of the simulated returns. Note that parameter combinations resulting
#' in non-stationary of the mean or variance will produce NAN's and that the maximum lag allowed for
#' any given parameter is 1.
#'
#' The default is ts_model = NULL, in which case time series time series properties are not induced, however if
#' ts_model = list() then the default values are list(omega = 5e-04, alpha = 0, gamma = NULL, beta = 0, mu = 0,
#' ar = NULL, ma = NULL, delta = 2). In order to set different parameters for each asset simply insert a vector
#' of length equal to the number of assets, the first element of the vector will correspond to Asset_1, the 2nd
#' to Asset_2 ect...
#'
#' See the sim_garch function's documentation and the "model" parameter under fGarch::garchSpec() for more details
#' regarding the parameters themselves.
#' @return a tidy tibble containing a date, Universe, Asset and Return column.
#'
#' @importFrom tidyr gather
#' @importFrom dplyr progress_estimated
#' @import dplyr
#' @import purrr
#'
#' @examples
#'
#' \dontrun{
#'
#' library(tidyverse)
#'
#' ### creating a correlation matrix to use as input in sim_asset_market
#' corr <- gen_corr(N = 20, Clusters = "none")
#'
#' ### simulating 550 periods of returns across 50 assets
#' set.seed(12542)
#' market_data <- sim_asset_market(corr,
#'                                 k = 550,
#'                                 mv_dist = "norm",
#'                                 left_cop_weight = 0.1,
#'                                 marginal_dist = "norm",
#'                                 ts_model = list(mu = 0.000002,
#'                                                 omega = 0.00005,
#'                                                 alpha = 0.098839,
#'                                                 beta = 0.899506,
#'                                                 ar = 0.063666,
#'                                                 ma = NULL,
#'                                                 gamma = 0.12194,
#'                                                 delta = 1.85))
#'
#'  ### Visualising the market
#'  market_data %>% group_by(Asset) %>%
#'  mutate(cum_ret = 100*cumprod(1 + Return)) %>%
#'          ggplot() +
#'          geom_line(aes(x = date, y = cum_ret, color = Asset)) +
#'          facet_wrap(~Asset) +
#'          theme(legend.position = "none")
#'
#' }
#' @export
mc_market <- function(corr,
                       N = 100,
                       k = 252,
                       mv_dist = "t",
                       mv_df = 3,
                       left_cop_weight = 0,
                       left_cop_param = 4,
                       marginal_dist = "norm",
                       marginal_dist_model = NULL,
                       ts_model = list()
) {

    1:N %>%
        map_dfr(~sim_market(corr = corr,
                            k = k,
                            mv_dist = mv_dist,
                            mv_df = mv_df,
                            left_cop_weight = left_cop_weight,
                            left_cop_param = left_cop_param,
                            marginal_dist = marginal_dist,
                            marginal_dist_model = marginal_dist_model,
                            ts_model = ts_model),
                .id = "Universe")
}
