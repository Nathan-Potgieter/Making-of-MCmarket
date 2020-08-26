README
================

# INTRO

This is the README file for Nathan Potgieter’s financial econometrics
project.

## Aim

The aim of this project is to develop a general and easy to use Monte
Carlo package that generates asset return data with a prespecified
correlation structure and dynamic dependencies. The user will be able to
adjust a “leverage” parameter, which determines the likelihood of
entering a “crisis period” characterized by extreme joint drawdowns. The
data will also be generated to exhibit ARIMA(p,q) + GARCH(q,p) features,
the parameters of which can be adjusted to induce alternative forms of
risk. Elliptical copulas are used to induce the correlation in the
simulated data, while Archmedian copulas are used to adjust the
likelihood of joint drawdowns. Various ARIMA(p,q) + GARCH(q,p)
structures can be called on to induce mean and variance persistence.

``` r
library(pacman)
p_load(tidyverse, copula, fGarch, libridate, bizdays)
```

    ## Installing package into 'C:/Users/Acer/Documents/R/win-library/4.0'
    ## (as 'lib' is unspecified)

    ## Warning: package 'libridate' is not available (for R version 4.0.2)

    ## Warning: unable to access index for repository http://www.stats.ox.ac.uk/pub/RWin/bin/windows/contrib/4.0:
    ##   cannot open URL 'http://www.stats.ox.ac.uk/pub/RWin/bin/windows/contrib/4.0/PACKAGES'

    ## Warning: 'BiocManager' not available.  Could not check Bioconductor.
    ## 
    ## Please use `install.packages('BiocManager')` and then retry.

    ## Warning in p_install(package, character.only = TRUE, ...):

    ## Warning in library(package, lib.loc = lib.loc, character.only = TRUE,
    ## logical.return = TRUE, : there is no package called 'libridate'

    ## Warning in p_load(tidyverse, copula, fGarch, libridate, bizdays): Failed to install/load:
    ## libridate

# Generating Covarience matrix

In this section I developed a simple function that allows the user to
easily generate a covarience matrix with the desired cluster structure.
Note that the majority of the code was writen by Nico Katzke. The
function is located in the gcVar.R code file.

### gcVar’s arguments

1.  N - is the number of assets in the universe

2.  Clusters - a character string specifying the type of cluster
    structure. Available options are “none”, for a correlation matrix
    with no clusters, “non-overlapping” for a correlation matrix with
    number one layer of clusters, and “overlapping” for a correlation
    matrix with Num\_Layers and Num\_clusters per layer.

3.  Num\_Clusters - if Clusters is equal to “non-overlapping” or “none”
    then Num\_Clusters is an integer value specifying the number of
    clusters. If Clusters = “overlapping” then Num\_Clusters must be a
    vector of length equal to Num\_Layers specifying the number of
    clusters per layer.

4.  Num\_Layers - an integer value between 1 and 4, specifying the
    number of cluster layers. Only needed of using “overlapping”
    clusters.

<!-- end list -->

``` r
#Co-Varience matrix generatimg function

gcVar <- 
  function(N = 50, Clusters = c("none", "non-overlapping", "overlapping") , Num_Clusters = NULL, Num_Layers = NULL){
    
Grps <- Num_Clusters
#set.seed(123)
    
if(Clusters == "none"){
    # Unclustered covariance matrix
    Sigma <- diag(N)
    for (i in 1:N) for (j in 1:N) Sigma[i,j] <- 0.9^abs(i-j)
    Sigma <- propagate::cor2cov(Sigma, runif(N, 1, 5))
    corr <- cov2cor(Sigma)
} else

if(Clusters == "non-overlapping"){
    #----------------------
    # distinct non-overlapping clusters:
    #----------------------
    
    if(is.null(Num_Clusters)) stop("Please provide a valid Num_Clusters argument when using Overlapping clusters")
    
    
    Sigma <- matrix(0.9, N, N)
    diag(Sigma) <- 1

    
for (i in 1:Grps) {
      ix <- seq((i-1) * N / Grps + 1, i * N / Grps)
      Sigma[ix, -ix] <- 0.0001                       #think about
    }
    Sigma <- propagate::cor2cov(Sigma, runif(N, 1, 5))
    corr <- cov2cor(Sigma)
} else
  
if(Clusters == "overlapping"){
    #----------------------
    # distinct overlapping clusters:
    #----------------------
  
  if(is.null(Num_Layers)|Num_Layers<2){
      stop("Please provide a valid Num_Layers argument when using Overlapping clusters")
      }else
  if(length(Num_Clusters) != Num_Layers){
      stop("Please provide a Num_Clusters argument with length equal to Num_Layers")
  }
    
  
    Sigma <- matrix(0.9, N, N)
    diag(Sigma) <- 1

    for (i in 1:Grps[1]) {
      ix <- seq((i-1) * N / Grps[1] + 1, i * N / Grps[1])
      Sigma[ix, -ix] <- 0.7
    }
    
    if(Num_Layers>=2){
        for (i in 1:Grps[2]) {
          ix <- seq((i-1) * N / Grps[2] + 1, i * N / Grps[2])
          Sigma[ix, -ix] <- 0.05
        } }else
    if(Num_Layers>=3){
        for (i in 1:Grps[3]) {
      ix <- seq((i-1) * N / Grps[3] + 1, i * N / Grps[3])
      Sigma[ix, -ix] <- 0.03
        } }else
    if(Num_Layers>=4){
        for (i in 1:Grps[4]) {
      ix <- seq((i-1) * N / Grps[4] + 1, i * N / Grps[4])
      Sigma[ix, -ix] <- 0
        } } 
    }

    Sigma <- propagate::cor2cov(Sigma, runif(N, 1, 5))  #Is this necessary???
    corr <- cov2cor(Sigma)

return(corr)

  }
```

Demonstrating the use of gcVar

## GeneratingRrandom Draws with numerous Copula Functions

### Elliptal copulas

Elliptal copulas such as the Gaussian and the student t copulas, allow
us to specify a correlation matrix as a parameter. Doing so allows one
to produce random draws of uniformly distributed variables, that contain
the correlation structure and joint distribution specified by the
copula. The chunk of code below demonstrates this functionality.

Unfortunately, both Elliptal copulas cannot be calibrated to exhibit
increased co-movements within the tails of the distribution. Therefore,
in the next section we examine some properties of Archimedean copulas.

### Archimedean copulas

Archimedean copulas such as the clayton, frank, gumbel and joe exhibit
increased dependence at the tails of the multivariate distribution. In
this section we will examine the clayton, …. copulas due to them
exhibiting enhanced left-tail dependencies. We will also have a look at
the hybrid BB1-BB6 which in which exibit increased dynamic dependenies
in both tails.

# Looking at some hybrid copulas

Tawn’s (1988) Theorem: Shows that a copula is a convex set and every
convex combination of existing copula functions is again a copula. See
“Extreme Dependence Structures and the Cross-Section of Expected Stock
Returns” page 8 & 9.

## hycop

The function below generates randomly distributed numbers from a hybrid
t and clayton copula. Need to think about how to calibrate df and
claycop parameters.

Arguments

  - Corr this is a correlation matrix uses as the parameter for the
    elliptical copula
  - elliptal\_copula family name of elliptal copula. Default is to use
    “t”, but “norm” is also accepted
  - df\_ellip a positive integer specifying the degrees of freedom for
    the student t elliptic copula. Only required when using
    elliptal\_copula = “t”.
  - left\_cop\_param a positive integer specifying the parameter of the
    **Clayton** copula.
  - left\_cop\_weight a value between 0 and 1 corresponding to the
    weight assigned to the left copula, when generating random draws
    from a hybrid copula.
  - marginal\_dist a character string specifying the marginal
    distribution of the simulated data. Must be “norm” or “t”, with the
    default generating uniformly distributed marginals.
  - df\_marginal\_dist a positive integer specifying the degrees of
    freedom parameter of the “t” distributed marginals.

<!-- end list -->

``` r
hycop <- function(corr,
                  elliptal_copula = c("norm", "t"), 
                  df_ellip = NULL, 
                  left_cop_param = 5,
                  left_cop_weight = 0.5,
                  T = 251, 
                  marginal_dist = NULL,
                  sd_marginal_dist = NULL,
                  nu_marginal_dist = NULL){
  
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
          
          if(is.null(sd_marginal_dist)|is.null(nu_marginal_dist))stop('Please supply a valid sd and nu parameter when using marginal_dist=="std". ')
          data <- apply(data, 2, qstd, sd = sd_marginal_dist, nu = nu_marginal_dist)
          
        }else
          if(marginal_dist=="norm"){
            
            data <- apply(data, 2, qnorm)
            
            }
              
              return(data)
            
}
```

Testing hycop

``` r
source("code/hycop.R")
source("code/gcVar.R")
set.seed(123)

Corr <- gcVar(N = 50, Clusters = "overlapping", Num_Layers = 3, Num_Clusters = c(10,5,2))
Corr %>% corrplot::corrplot()
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
#Using marginal_dist="std"
set.seed(123)
data <- hycop(Corr, elliptal_copula = "t", df_ellip = 10, 
              left_cop_param = 6, left_cop_weight = 0.5, T = 10000,
              marginal_dist = "std", nu_marginal_dist = 8, sd_marginal_dist = 1)

#Notice how the distribution change when moving to variables outside the specified risk cluster
data %>% plot(main = 'Hycop: std-distmarginals, nu = 8,  sd = 1') 
```

![](README_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

``` r
plot(data[,1], data[,10], main = 'Hycop: std-distmarginals, nu = 8,  sd = 1')
```

![](README_files/figure-gfm/unnamed-chunk-2-3.png)<!-- -->

``` r
plot(data[,1], data[,ncol(data)], main = 'Hycop: t-dist marginals, nu = 8, sd = 1')
```

![](README_files/figure-gfm/unnamed-chunk-2-4.png)<!-- -->

``` r
#Note how the correlation matrix has been maintained
data %>% cor %>%  corrplot::corrplot()
```

![](README_files/figure-gfm/unnamed-chunk-2-5.png)<!-- -->

# Looking at options for marginal distributions

The std distribution.

``` r
# Display the Student's t distributions with various
# degrees of freedom and compare to the normal distribution
std.dist <- function(xlim = c(-4,4), ylim = c(0,0.7), param = c(3, 5, 7, 30), legend = "topright"){
  library(fGarch)
x <- seq(xlim[1], xlim[2], length=100)
hx <- dnorm(x)

degf <- param
colors <- c("steelblue", "blue", "darkgreen", "gold", "black")
labels <- c(glue::glue("nu={degf}"), "normal")

plot(x, hx, type="l", lty=2, xlab="x value",
     ylab="Density", main="Comparison of std Distributions", ylim=ylim)

for (i in 1:length(param)){
  lines(x, dstd(x,nu=degf[i], mean = 0, sd = 1), lwd=2, col=colors[i])
}

legend(legend, inset=.05, title="Distributions",
       labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)
}

#std vs normal
std.dist(param = seq(3, 9, length.out = 4), ylim = c(0,0.7))
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
#std vs normal left tail
std.dist(param = seq(3, 9, length.out = 4), ylim = c(0,0.0025), xlim = c(-10,-4), legend = "topleft")
```

![](README_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

Looking at the skewed generalizd t distribution

``` r
p_load(sgt)
#-----------------------------------------
#Looking at the Skewnes parameter in sgt distribution
#-----------------------------------------

x <- seq(-4, 4, length=100)
hx <- dnorm(x)

degf <- -c(seq(0,0.45,length = 4))
colors <- c("red", "blue", "darkgreen", "gold", "black")
labels <- c(glue::glue("Param = {degf}"), "normal")

plot(x, hx, type="l", lty=2, xlab="x value",
  ylab="Density", main="Comparison of SGD Distributions with different lambda's", ylim = c(0,0.5))

for (i in 1:4){
  lines(x, dsgt(x,lambda = degf[i]), lwd=2, col=colors[i])
}

legend("topright", inset=.05, title="Distributions",
  labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
#-----------------------------------------
#now looking at p and q parameters (kurtosis) with lambda = 0.2
#-----------------------------------------

#altering p
x <- seq(-4, 4, length=100)
hx <- dnorm(x)

degf <- c(seq(1.01,2,length = 4))
colors <- c("red", "blue", "darkgreen", "gold", "black")
labels <- c(glue::glue("Param = {degf}"), "normal")

plot(x, hx, type="l", lty=2, xlab="x value",
  ylab="Density", main="Comparison of SGD Distributions with different p parameters", ylim = c(0,0.8))

for (i in 1:4){
  lines(x, dsgt(x,p = degf[i]), lwd=2, col=colors[i])
}

legend("topright", inset=.05, title="Distributions",
  labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)
```

![](README_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

``` r
#altering q
x <- seq(-4, 4, length=100)
hx <- dnorm(x)

degf <- c(seq(1.5,4,length = 4)) %>% round(digits = 2)
colors <- c("red", "blue", "darkgreen", "gold", "black")
labels <- c(glue::glue("Param = {degf}"), "normal")

plot(x, hx, type="l", lty=2, xlab="x value",
  ylab="Density", main="Comparison of SGD Distributions with different p parameters", ylim = c(0,0.8))

for (i in 1:4){
  lines(x, dsgt(x,q = degf[i]), lwd=2, col=colors[i])
}

legend("topleft", inset=.05, title="Distributions",
  labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)
```

![](README_files/figure-gfm/unnamed-chunk-4-3.png)<!-- -->

# Introducing autocorrelation and Volitility clustering

In this step I introduce autocorrelation and volatility using an
AR(p,q)+GARCH(q,p) model.

  - “The leptokurtosis, clustering volatility and leverage effects
    characteristics of financial time series justifies the GARCH
    modeling approach. The non-linear characteristic of the time series
    is used to check the Brownian motion and investigate into the
    temporal evolutionary patterns. The nonlinear methods of forecasting
    and signal analysis are gaining popularity in stock market because
    of their robustness in feature extraction and classiﬁcation.”
    <https://towardsdatascience.com/garch-processes-monte-carlo-simulations-for-analytical-forecast-27edf77b2787>

My go at writing a GARCHSIM function, much of the code was borrowed from
fGarch’s garchspec and garchsim functions.

``` r
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
    colnames(data) <- glue::glue("Asset:{1:ncol(data)}")
    data <- data[-(1:m), ]   #removes burn in data
    
    
    if(simple == TRUE){
      data[,3]
      }else data
  
}
```

## Simulating an asset market

``` r
source("code/gcVar.R")
source("code/hycop.R")
source("code/Garch.sim.R")

set.seed(1234)

corr <- gcVar(N = 20, Clusters = "none", Num_Layers = 3, Num_Clusters = c(2,4,5))  #need to work on Corr matrix generation


inno <- hycop(corr, elliptal_copula = "t",
                      df_ellip = 4, 
                      left_cop_param = 5, 
                      left_cop_weight =0,
                      T = 500,
                      marginal_dist = "norm", 
                      sd_marginal_dist = 1,
                      nu_marginal_dist = 3)
colnames(inno) <- glue::glue("Asset_{1:ncol(inno)}") #move to hycop
#inno  %>% cor() %>% corrplot::corrplot()

model <- list(mu = 0.000048,        #Parameters from Statistics and Data Analysis for Financial Engineering pg.421-423
              omega = 0.000050, 
              alpha = 0.098839, 
              beta = 0.899506, 
              ar = 0.063666,
              ma = NULL,
              gamma = 0.12194,
              delta = 1.476643)

#Making Simdat tidy
simdat <- suppressMessages(
  1:ncol(inno) %>% map(~Garch.sim(model, rev(inno[,.x]))) %>% reduce(bind_cols)
  )
colnames(simdat) <- glue::glue('Asset_{1:ncol(simdat)}')

tidy_simdat <- simdat %>%  mutate(date = 1:nrow(simdat)) %>% 
  gather(key = Asset, value = Return, -date)

#Correlation matrix of GARCH sim data
#simdat %>% spread(key = Asset, value = Return) %>% select(-date) %>% cor() %>% corrplot::corrplot()

#generating and plotting Cum returns
tidy_simdat %>%  
  arrange(date) %>% 
  group_by(Asset) %>%
  mutate(Cum_Return = cumprod(1 + Return)) %>%
  ggplot() + geom_line(aes(x=date, y=Cum_Return, color = Asset)) +
  facet_wrap(~Asset, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "none") 
```

![](README_files/figure-gfm/sim%20market-1.png)<!-- -->

``` r
tidy_simdat %>% ggplot() +
  geom_line(aes(x=date,y=Return, color = Asset)) +
  facet_wrap(~Asset) +
  theme_bw() +
  theme(legend.position = "none")
```

![](README_files/figure-gfm/sim%20market-2.png)<!-- -->

### Testing if garchsim has the same issue where vol is exceeding large at beginning of simulations

``` r
set.seed(1234)
spec <- garchSpec(model = list(mu = 0.000048,
                               omega = 0.000050, 
                               alpha = 0.098839, 
                               beta = 0.899506, 
                               ar = 0.063666,
                               ma = NULL,
                               gamma = 0.12194,
                               delta = 1.476643))


garchSim(spec, n = 500) %>% plot()
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

## Simulating an ensemble of asset markets: Prototype

``` r
Sim.Asset.Market <- function(corr, T = 500){    #note that length out is 499 not 500,  look at garchsim to fix
  
  inno <- hycop(corr, elliptal_copula = "t",
                      df_ellip = 4, 
                      left_cop_param = 5, 
                      left_cop_weight = 0,
                      T = 500,
                      marginal_dist = "std", 
                      sd_marginal_dist = 1,
                      nu_marginal_dist = 9)

model <- list(mu = 0.002,
              omega = 1.76E-04, 
              alpha = 0.053, 
              beta = 0.922, 
              ar = 0.01,
              ma = NULL)

#Applying garch.sim to each column in simdat
simdat <- suppressMessages(
  1:ncol(inno) %>% map(~Garch.sim(model, inno[,.x])) %>% reduce(bind_cols) )
  
# Naming Each coloumn
colnames(simdat) <- glue::glue('Asset_{1:ncol(simdat)}')

#adding a date column
bizdays::create.calendar(name = "weekdays", weekdays = c("saturday", "sunday"))
start_date <- Sys.Date()
end_date <- bizdays::offset(start_date, nrow(simdat), "weekdays")
weekdays <- rmsfuns::dateconverter(start_date, end_date, "weekdays")[1:nrow(simdat)]


simdat <- simdat %>% mutate(date = weekdays, .before = `Asset_1`)
simdat %>% gather(key=Asset, value = Return, -date)
}
```

Now lets use the function above to run our first MC simulation.

``` r
#Is there a better way to name cols? i.e before we reduce?
mc.data <- suppressMessages(
1:10 %>% map(~Sim.Asset.Market(corr, T = 500)) %>% reduce(left_join, by = c("date","Asset"))
)
colnames(mc.data) <- c("date", "Asset", glue::glue("Universe_{1:(ncol(mc.data)-2)}"))

#This is how I want my final output to look. 
mc.data %>% gather(Universe, Return, c(-date,-Asset))
```

    ## # A tibble: 99,800 x 4
    ##    date       Asset   Universe     Return
    ##    <date>     <chr>   <chr>         <dbl>
    ##  1 2020-08-25 Asset_1 Universe_1  0.0187 
    ##  2 2020-08-26 Asset_1 Universe_1 -0.0298 
    ##  3 2020-08-27 Asset_1 Universe_1 -0.0716 
    ##  4 2020-08-28 Asset_1 Universe_1 -0.0436 
    ##  5 2020-08-31 Asset_1 Universe_1 -0.00279
    ##  6 2020-09-01 Asset_1 Universe_1 -0.0180 
    ##  7 2020-09-02 Asset_1 Universe_1  0.0417 
    ##  8 2020-09-03 Asset_1 Universe_1  0.120  
    ##  9 2020-09-04 Asset_1 Universe_1  0.0960 
    ## 10 2020-09-07 Asset_1 Universe_1  0.0142 
    ## # ... with 99,790 more rows
