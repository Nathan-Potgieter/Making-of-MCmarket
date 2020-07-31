README
================

# INTRO

This is the README file for Nathan Potgieter’s financial econometrics
project.

The aim of this project is to develop a general and easy to use Monte
Carlo package, that generates financial data with the possibility of
extreme joint down movements, as observed during financial crisis.

``` r
library(tidyverse)
```

    ## -- Attaching packages ------------------------------------------------------------------------------------------------------------------------ tidyverse 1.3.0 --

    ## v ggplot2 3.3.2     v purrr   0.3.4
    ## v tibble  3.0.3     v dplyr   1.0.0
    ## v tidyr   1.1.0     v stringr 1.4.0
    ## v readr   1.3.1     v forcats 0.5.0

    ## -- Conflicts --------------------------------------------------------------------------------------------------------------------------- tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

## Generating Covarience matrix

In this section I developed a simple function that allows the user to
easily generate a covarience matrix with the desired cluster structure.
Note that the majority of the code was writen by Nico Katzke. The
function is located in the gcVar.R code file.

##### Question refering to this chunk of code

1.  How to I specify the function argument so that; if(Clusters ==
    “overlapping”) then Num,Clusters must be a vector of length 3.
    that way the user can specify the number of clusters at each layer.

2.  
<!-- end list -->

``` r
#Co-Varience matrix generatimg function

gcVar <- function(N = 50, Clusters = c("none", "non-overlapping", "overlapping") , Num.Clusters = 10){

#If Clusters = "overlapping" then Num.Clusters must be a vector fo length 3
    
N <- N
Grps <- Num.Clusters
#set.seed(123)
    
if(Clusters == "none"){
    # Unclustered covariance matrix
    Sigma <- diag(N)
    for (i in 1:N) for (j in 1:N) Sigma[i,j] <- 0.9^abs(i-j)
    Sigma <- propagate::cor2cov(Sigma, runif(N, 1, 5))
    corr <- cov2cor(Sigma)
}

if(Clusters == "non-overlapping"){
    #----------------------
    # distinct non-overlapping clusters:
    #----------------------
    Sigma <- matrix(0.9, N, N)
    diag(Sigma) <- 1
    
    for (i in 1:Grps) {
      ix <- seq((i-1) * N / Grps + 1, i * N / Grps)
      Sigma[ix, -ix] <- 0.1
    }
    Sigma <- propagate::cor2cov(Sigma, runif(N, 1, 5))
    corr <- cov2cor(Sigma)
}

if(Clusters == "overlapping"){
    #----------------------
    # distinct overlapping clusters:
    #----------------------
    Sigma <- matrix(0.9, N, N)
    diag(Sigma) <- 1
    
    Grps <- 10  #Grps[1]
    for (i in 1:Grps) {
      ix <- seq((i-1) * N / Grps + 1, i * N / Grps)
      Sigma[ix, -ix] <- 0.6
    }
    
    Grps <- 5   #Grps[2]
    for (i in 1:Grps) {
      ix <- seq((i-1) * N / Grps + 1, i * N / Grps)
      Sigma[ix, -ix] <- 0.3
    }
    
    Grps <- 2   #Grps[3]
    for (i in 1:Grps) {
      ix <- seq((i-1) * N / Grps + 1, i * N / Grps)
      Sigma[ix, -ix] <- 0.1
    }
    
    Sigma <- propagate::cor2cov(Sigma, runif(N, 1, 5))
    corr <- cov2cor(Sigma)
}

return(corr)

}
```

``` r
gcVar(N = 50, Clusters = "none") %>% corrplot::corrplot()
```

![](README_files/figure-gfm/using%20gcVar-1.png)<!-- -->

``` r
gcVar(N = 50, Clusters = "non-overlapping") %>% corrplot::corrplot()
```

![](README_files/figure-gfm/using%20gcVar-2.png)<!-- -->

``` r
gcVar(N = 50, Clusters = "overlapping") %>% corrplot::corrplot()
```

![](README_files/figure-gfm/using%20gcVar-3.png)<!-- -->

## GeneratingRrandom Draws with numerous Copula Functions
