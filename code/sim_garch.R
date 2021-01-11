#' @title sim_garch
#' @description This function takes a vector of random numbers, referred to as innovations, and
#' induces mean and variance persistence by inserting them into an ARMA(1,1) + APARCH(1,1) model.
#' @param innovations a vector containing the random numbers/ the innovations of the
#' ARIMA + GARCH process.
#' @param omega a positive value defining the coefficient of the variance equation, default is 5e-04.
#' @param gamma a value defining the APARCH leverage parameter in the variance equation. The default
#' of 0, implies no leverage effect and therefore corresponds with the standard GARCH model.
#' @param alpha a value defining the value of the autoregressive variance coefficient, default is 0.
#' @param beta a value defining the variance coefficient, default is 0.
#' @param mu  a value defining the mean, default is 0.
#' @param ar  a value defining the autoregressive ARMA coefficient, default is 0.
#' @param ma a value defining the moving average ARMA coefficient, default is 0.
#' @param delta a strictly positive value the delta parameter of the APARCH model. The default is 2,
#' which corresponds with the standard GARCH model.
#' @param simple a logical parameter indicating if the output should be a simple vector containing just the
#' resulting ARIMA + GARCH series, or if FALSE a three column dataframe containing z - the innovations, h - the
#'  conditional variance and y - ARMA + APARCH series.
#' @note  (1) It is suggested that the randomly distributed numbers be mean zero and standard
#' deviation one, as these attributes can be set by the model argument.
#'
#' (2) For more information on the ARMA + APARCH parameters see:
#'
#' Ruppert, D. and Matteson, D.S., 2011. Statistics and data analysis for financial engineering (Vol. 13). New York: Springer.
#'
#'  @return if simple = TRUE a vector of the resulting ARMA + APARCH series, else if simple = FALSE a
#' three column dataframe containing z - the innovations, h - the conditional variance and y - ARMA +
#' APARCH series. Note the length of the resulting series will one observation less than that that of the innovations
#' as ARMA(1,1) + APARCH(1,1) model effectively consumes this lag when producing its first value.
#'
#' @importFrom dplyr tibble
#'
#' @examples
#'
#' \dontrun{
#'
#' library(tidyverse)
#'
#' ### creating series of 501 innovations
#' inno <-  rnorm(501)
#'
#' ### This produces a ARMA + APARCH series of length 500.
#' GARCH <- sim_garch(inno,
#'                    mu = 0.000002,
#'                    omega = 0.00005,
#'                    alpha = 0.098839,
#'                    beta = 0.899506,
#'                    ar = 0.063666,
#'                    ma = NULL,
#'                    gamma = 0.12194,
#'                    delta = 1.85,
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

sim_garch <- function(innovations,
                      omega = NULL, alpha = NULL, gamma = NULL, beta = NULL,
                      mu = NULL, ar = NULL, ma = NULL, delta = NULL,
                      simple = TRUE) {

    # Creating model list with defaults
    model <- list(omega = ifelse(is.null(omega), 5e-04, omega),
                  alpha = ifelse(is.null(alpha), 0, alpha),
                  gamma = ifelse(is.null(gamma), 0, gamma),
                  beta = ifelse(is.null(beta), 0, beta),
                  mu = ifelse(is.null(mu), 0, mu),   #changed form NULL to 0
                  ar = ifelse(is.null(ar), 0, ar),
                  ma = ifelse(is.null(ma), 0, ma),
                  delta = ifelse(is.null(delta), 2, delta))


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
    n <- length(innovations) - 1

    #Generating innovations
    z_length <- n + max.order
    z <- c(innovations)[(2-max.order):length(innovations)]  #z is the vector of random innovation

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
        data <- c(NA, y[(m + 1):(n + m)]) #removes burn in data and only provides ARIMA + GARCH series.
    } else {
        data <- tibble(z = c(NA, z[(m + 1):(n + m)]),  # Innovations
                       h = c(NA, h[(m + 1):(n + m)]^deltainv), #Conditional Standard deviation
                       y = c(NA, y[(m + 1):(n + m)]))  # ARIMA + GARCH series
    }
    return(data)
}
