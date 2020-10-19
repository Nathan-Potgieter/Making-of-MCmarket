#Co-Varience matrix generatimg function

gen_corr <-
    function(N = 50, Clusters = c("none", "non-overlapping", "overlapping") , Num_Clusters = NULL, Num_Layers = NULL){

        Grps <- Num_Clusters
        #set.seed(123)

        if(Clusters == "none"){
            # Unclustered covariance matrix
            Sigma <- diag(N)
            for (i in 1:N) for (j in 1:N) Sigma[i,j] <- 0.7^abs(i-j)
            Sigma <- propagate::cor2cov(Sigma, runif(N, 1, 5))
            corr <- cov2cor(Sigma)
        } else

            if(Clusters == "non-overlapping"){
                #----------------------
                # distinct non-overlapping clusters:
                #----------------------

                if(is.null(Num_Clusters)) stop("Please provide a valid Num_Clusters argument when using non-overlapping clusters")


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