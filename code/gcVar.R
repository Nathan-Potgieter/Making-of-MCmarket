#Co-Varience matrix generatimg function

gcVar <- function(N = 25, Runs = 251, Clusters = 10, hivol.weight = runif(N, 1, 5)){

    #----------------------
    # distinct overlapping clusters:
    #----------------------

    N <- N
    Runs <- Runs


    #thetas <- seq(0, 1, by = 0.05)
    #thetas <- seq(1,50,by=2)
    #thetas[length(thetas)] <- 50
    #thetas <- (thetas - 1)/49

    # Unclustered covariance matrix
    Sigma <- diag(N)
    for (i in 1:N) for (j in 1:N) Sigma[i,j] <- 0.9^abs(i-j)
    Sigma <- matrix(0.9, N, N)
    diag(Sigma) <- 1

    Grps <- Clusters
    for (i in 1:Grps) {
        ix <- seq((i-1) * N / Grps + 1, i * N / Grps)
        Sigma[ix, -ix] <- 0.6
    }
    Sigma <- propagate::cor2cov(Sigma, Var)
}