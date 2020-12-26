#' @title gen_corr
#' @description This function generates ad hoc correlation matrices with a set number of clusters and up to 4 layers.
#' @param D The number of variables, generates an N by N correlation matrix.
#' @param Clusters a character string specifying the type of cluster structure. Available options are "none", for a correlation matrix with no clusters, "non-overlapping" for a correlation matrix with one layer of clusters, and "overlapping" for a correlation matrix with up to 4 layers and a set number of clusters per layer.
#' @param Num_Clusters if Clusters = "non-overlapping" or Clusters = "none" then Num_Clusters is an integer value specifying the number of clusters. If Clusters = "overlapping" then Num_Clusters must be a vector of length equal to Num_Layers specifying the number of clusters per layer.
#' @param Num_Layers an positive integer value between 1 and 4, specifying the number of cluster layers. Only needed if using "overlapping" clusters.
#' @return a N by N correlation matrix.
#'
#' @examples
#'
#' \dontrun{
#' library(ggcorrplot)
#'
#' ### This generates a 50 by 50 correlation matrix with no clusters.
#' gen_corr(N = 50, Clusters = "none) %>%
#'          ggcorrplot(title = "Overlapping Clusters")
#'
#' ### This generates a 50 by 50 correlation matrix with 5 non-overlapping clusters.
#' gen_corr(N = 50, Clusters = "non-overlapping) %>%
#'          ggcorrplot(title = "Overlapping Clusters")
#'
#' ### This generates a 60 by 60 correlation matrix consisting
#' ### of 4 layers with 10, 5, 3 and 2 clusters respectively.
#' gen_corr(N = 60,
#'          Clusters = "overlapping",
#'          Num_Layers = 4,
#'          Num_Clusters = c(10,5,3,2)) %>%
#'                   ggcorrplot::ggcorrplot(title = "Overlapping Clusters")
#'
#' }
#' @export
gen_corr <- function (D = 50,
                      Clusters = c("none", "non-overlapping", "overlapping"),
                      Num_Clusters = NULL,
                      Num_Layers = NULL) {

        Grps <- Num_Clusters
        #set.seed(123)

        if(Clusters == "none"){
                # Unclustered covariance matrix
                Sigma <- diag(D)
                for (i in 1:D) for (j in 1:D) Sigma[i,j] <- 0.9^abs(i-j)
                Sigma <- propagate::cor2cov(Sigma, runif(D, 1, 5))
                corr <- cov2cor(Sigma)
        } else

                if(Clusters == "non-overlapping"){
                        #-----------------------------------
                        # distinct non-overlapping clusters:
                        #-----------------------------------

                        if(is.null(Num_Clusters)) stop("Please provide a valid Num_Clusters argument when using Overlapping clusters")


                        Sigma <- matrix(0.6, D, D)
                        diag(Sigma) <- 1


                        for (i in 1:Grps) {
                                ix <- seq((i-1) * D / Grps + 1, i * D / Grps)
                                Sigma[ix, -ix] <- 0.0001                       #think about
                        }
                        Sigma <- propagate::cor2cov(Sigma, runif(D, 1, 5))
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


                                Sigma <- matrix(0.8, D, D)
                                diag(Sigma) <- 1

                                for (i in 1:Grps[1]) {
                                        ix <- seq((i-1) * D / Grps[1] + 1, i * D / Grps[1])
                                        Sigma[ix, -ix] <- 0.5
                                }
                                if(Num_Layers>=2){
                                        for (i in 1:Grps[2]) {
                                                ix <- seq((i-1) * D / Grps[2] + 1, i * D / Grps[2])
                                                Sigma[ix, -ix] <- 0.3
                                        } }
                                if(Num_Layers>=3){
                                        for (i in 1:Grps[3]) {
                                                ix <- seq((i-1) * D / Grps[3] + 1, i * D / Grps[3])
                                                Sigma[ix, -ix] <- 0.15
                                        } }
                                if(Num_Layers>=4){
                                        for (i in 1:Grps[4]) {
                                                ix <- seq((i-1) * D / Grps[4] + 1, i * D / Grps[4])
                                                Sigma[ix, -ix] <- 0.05
                                        } }
                        }

        Sigma <- propagate::cor2cov(Sigma, runif(D, 1, 5))
        corr <- cov2cor(Sigma)

        return(corr)

}
