#' @title sim_garch
#' @description This function takes a vector of random numbers and induces mean and variance persistence by plugging them into an ARIMA + GARCH model as the innovations.
#' @note  It is suggested that the randomly distributed numbers be mean zero and standard deviation one, as those attributes are better set in the model argument.
#' @param innovations a vector containing the random numbers/ the innovations of the ARIMA + GARCH process.
#' @param model a list containing various ARIMA + GARCH, AP-GARCH, GJR-GARCH, ect... parameters allowing once to specify the time series properties of the simulated returns. Note that parameter combinations resulting in non-stationary of the mean or variance will produce NAN's. The default values are set as list(omega = 5e-04, alpha = 0, gamma = NULL, beta = 0, mu = 0, ar = NULL, ma = NULL, delta = 2). Note that omega is a key input indicating the constant coefficient of the variance equation and that each parameter can be a vector of length <= 5, this will induce a higher order ARIMA + GARCH process. For more information on the various parameters see the user manual.
#' @param simple a logical parameters if the output should be a simple vector containing the resulting ARIMA + GARCH series, or if FALSE a three column dataframe containing z - the innovations, h - the conditional variance and y - ARIMA + GARCH series.
#' @return if simple = TRUE a vector of the resulting ARIMA + GARCH series, else if simple = FALSE a a three column dataframe containing z - the innovations, h - the conditional variance and y - ARIMA + GARCH series. Note the length of the reasulting series will be less that that of the innovations as ARIMA + GARCH models require up to 5 innovations (depending on the max order) to produce the first result.
#'
#' @importFrom xts as.xts timeBased
#'
#' @import base ????????????????????????????????????????????????????????????????
#'
#' @importFrom rlang :=
#' @importFrom rlang := enquo quo_get_expr
#' @examples
#'
#' \dontrun{
#'
#' library(tidyverse)
#'
#' ### creating series of 501 innovations
#' inno <-  rnorm(501)
#'
#' ### This produces a ARIMA + GARCH series of length 500.
#' GARCH <- sim_garch(inno,
#'                    model = list(mu = 0.000002,
#'                                 omega = 0.00005,
#'                                 alpha = 0.098839,
#'                                 beta = 0.899506,
#'                                 ar = 0.063666,
#'                                 ma = NULL,
#'                                 gamma = 0.12194,
#'                                 delta = 1.85),
#'                    simple = FALSE)
#'
#'  ### Visualising the resulting series
#'  GARCH %>% mutate(period = 1:n()) %>%
#'        gather(key, value, -period) %>%
#'        ggplot() +
#'        geom_line(aes(x = period, y = value, color = key)) +
#'        facet_wrap(~key, scales = "free_y", ncol = 1)
#'
#' }
#' @export

sim_garch <- function(innovations, model= list(), simple = TRUE) {

    #default parameters for garch model
    default <- list(omega = 5e-04,
                    alpha = 0,
                    gamma = NULL,
                    beta = 0,
                    mu = 0,   #changed form NULL to 0
                    ar = NULL,
                    ma = NULL,
                    delta = 2)

    default[names(model)] <- model
    model <- default

    #obtaining parameters and lag orders from model object
    mu <- model$mu
    ar <- model$ar
    ma <- model$ma
    omega <-  model$omega
    alpha <- model$alpha
    gamma <- model$gamma
    beta <- model$beta
    delta <- model$delta
    deltainv <- 1/delta
    order.ar <- length(ar)
    order.ma <- length(ma)
    order.alpha <- length(alpha)
    order.beta <- length(beta)
    max.order <- max(order.ar, order.ma, order.alpha, order.beta)
    n <- length(innovations) - 5

    if(max.order>5) stop("Please supply a volitility model with max order less than or equal to 5")

    #Generating innovations
    z_length <- n + max.order
    z <- c(innovations)[1:z_length]  #z is the random innovation

    h <- c(rep(model$omega/(1 - sum(model$alpha) - sum(model$beta)),
               times = max.order), rep(NA, n))    #h is the conditional standard deviation

    y <- c(rep(model$mu/(1 - sum(model$ar)), times = max.order), rep(NA, n))  #Observations from the resulting garch process
    m <- max.order

    #Inducing ARIMA + GARCH
    eps <- h^deltainv * z  #this part often breaks depending on GARCH parameters chosen (must be stationary)

    for (i in (m + 1):(n + m)) {
        h[i] = omega + sum(alpha * (abs(eps[i - (1:order.alpha)]) -
                                        gamma * (eps[i - (1:order.alpha)]))^delta) +
            sum(beta * h[i - (1:order.beta)])

        eps[i] = h[i]^deltainv * z[i]
        y[i] = mu + sum(ar * y[i - (1:order.ar)]) + sum(ma * eps[i - (1:order.ma)]) + eps[i]
    }

    if(simple == TRUE) {
        data <- y[(m + 1):(n + m)] #removes burn in data and only provides ARIMA + GARCH series.
    } else {
        data <- tibble(z = z[(m + 1):(n + m)],  # Innovations
                       h = h[(m + 1):(n + m)]^deltainv, #Conditional Standard deviation
                       y = y[(m + 1):(n + m)])  # ARIMA + GARCH series
    }
    return(data)
}