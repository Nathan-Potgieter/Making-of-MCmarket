
Garch.sim <- function(model= list(), innovations, simple = TRUE){

    #default parameters for garch model
    default <- list(omega = 1e-06,
                    alpha = 0.1,
                    gamma = NULL,
                    beta = 0.8,
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
    n <- length(innovations)

    #Generating innovations
    z <- c(rnorm(max.order), innovations)  #must change this later

    h <- c(rep(model$omega/(1 - sum(model$alpha) - sum(model$beta)),
               times = max.order), rep(NA, n))    #sd's

    y <- c(rep(model$mu/(1 - sum(model$ar)), times = max.order), rep(NA, n))  #garch simulations
    m <- max.order

    #simulating GARCH
    eps <- h^deltainv * z  #this part often breaks depending on GARCH parameters chosen

    for (i in (m + 1):(n + m)) {
        h[i] = omega + sum(alpha * (abs(eps[i - (1:order.alpha)]) -
                                        gamma * (eps[i - (1:order.alpha)]))^delta) +
            sum(beta * h[i - (1:order.beta)])

        eps[i] = h[i]^deltainv * z[i]
        y[i] = mu + sum(ar * y[i - (1:order.ar)]) + sum(ma *
                                                            eps[i - (1:order.ma)]) + eps[i]
    }

    data <- cbind(z = z[(m + 1):(n + m)],
                  sigma = h[(m + 1):(n + m)]^deltainv,
                  y = y[(m + 1):(n + m)])

    rownames(data) <- as.character(1:n)  # may want to change to dates
    data <- data[-(1:m), ]   #removes burn in data

    if(simple == TRUE){
        data[,3]
    }else data

}
