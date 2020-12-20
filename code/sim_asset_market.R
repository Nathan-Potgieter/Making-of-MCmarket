#' @title sim_asset_market
#' @description This function produces a series of returns for an asset market with the given correlation matrix. The user can adjust the markets left tail dependency as well as the markets distribution and the univariate distributions of the returns.
#' @note  It is suggested that the marginal distributions be set to mean zero and standard deviation one. Those attributes are better set in the ts_model argument.
#' @param corr a correlation matrix that the simulated date will adhere to. Note that the number of variables simulated is equal to the number of columns in the correlation matrix.
#' @param k a positive integer indicating the number of time periods to simulate. Note that the number of periods generated is actually equal to k + 5 as these extra observations are needed when applying time series properties to the data.
#' @param mv_dist a string specifying the multivariate distribution. Can be one of c("norm", "t") referring to the multivariate normal and t distributions respectively. Default is 3.
#' @param mv_df degrees of freedom of the multivariate distribution, required when mv_dist = "t".
#' @param left_cop_weight a positive value between zero and one indicating the weight applied to the Clayton copula when creating the multivariate distribution. Note that a value between zero and one essentially generates a hybrid distribution between the chosen mv_dist and the Clayton copula. Therefore, the greater the left_cop_weight the less the data will reflect the correlation structure. Default is set to 0.
#' @param left_cop_param a positive value indicating the parameter of the Clayton copula. Default is 4.
#' @param marginal_dist a string variable specifying the univariate distribution of each variable. Can be one of c("norm", "t", "sgt") referring to the normal, student-t and skewed-generalized-t distributions respectively. Default is "norm".
#' @param  marginal_dist_model list containing the relevant parameters for the chosen marginal_dist. marginal_dist = "norm" accepts a mean and standard deviation with defaults list(mu = 0, sigma = 1) respectively. marginal_dist = "t" accepts the non-centrality and degrees of freedom arguments, default values are list(mu = 0, df = 5). marginal_dist = "sgt" accepts the mean, sd, lambda, p and q parameters list(mu = 0, sigma = 1, lambda, p, q). Note lambda, p and q have no defaults and must therefore be set by the user.
#' @param ts_model a list containing various ARIMA + APGARCH parameters allowing one to specify the time series properties of the simulated returns. Note that parameter combinations resulting in non-stationary of the mean or variance will produce NAN's and that the maximum lag allowed for any given parameter is 5. The default values are set as list(omega = 5e-04, alpha = 0, gamma = NULL, beta = 0, mu = 0, ar = NULL, ma = NULL, delta = 2). Note that omega is a key input indicating the constant coefficient of the variance equation and that each parameter can be a vector of length <= 5, this will induce a higher order ARIMA + GARCH process. For more information on the various parameters see the user manual.
#' @return a tidy tibble containing a date, Asset and Return column.
#'
#' @importFrom xts as.xts timeBased
#'
#' @import copula
#' @import glue
#' @import dplyr
#' @import purrr
#' @import sgt
#' @import rmsfuns
#'
#' @importFrom rlang :=
#' @importFrom rlang := enquo quo_get_expr
#' @examples
#'
#' \dontrun{
#'
#' library(tidyverse)
#'
#' ### creating a correlation matrix of 50 assets to use as an input in sim_asset_market.
#' corr <- gen_corr(N = 20, Clusters = "none")
#'
#' ### simulating 500 periods of returns across 50 assets.
#' set.seed(12345)
#' market_data <- sim_asset_market(corr,
#'                                 k = 500,
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

sim_asset_market <- function(corr,
                             k = 252,
                             mv_dist = "t",
                             mv_df = 3,
                             left_cop_weight = 0,
                             left_cop_param = 4,
                             marginal_dist = "norm",
                             marginal_dist_model = NULL,         # may want to change to a list
                             ts_model = list()
) {
    #Simulating innovations
    inno <- sim_inno(corr = corr,
                     mv_dist = mv_dist,
                     mv_df = mv_df,
                     left_cop_param = left_cop_param,
                     left_cop_weight = left_cop_weight,
                     marginal_dist = marginal_dist,
                     marginal_dist_model = marginal_dist_model,
                     k = k)

    #creating a date vector
    start_date <- Sys.Date()
    dates <- rmsfuns::dateconverter(StartDate = start_date,
                                    EndDate = start_date %m+% lubridate::days(k-1),
                                    Transform = "alldays")

    if (is.null(ts_model)) {
        return(
            inno[6:nrow(inno),] %>% mutate(date = dates, .before = `Asset_1`) %>%
                gather(key = Asset, value = Return, -date)
        )
    } else

    #Applying sim_garch to each column in simdat

    simdat <- inno %>% map_dfc(~sim_garch(innovations = .x, model =  ts_model))

    #Creating final df
    simdat %>%
        mutate(date = dates, .before = `Asset_1`) %>%
        gather(key = Asset, value = Return, -date)
}