#' @title gen_corr
#' @description This function allows users to easily generate ad hoc correlation
#' matrices with a set number of clusters and up to 4 layers.
#' @param D The number of variables, generates an D by D correlation matrix.
#' @param clusters a character string specifying the type of cluster structure.
#' Available options are "none", for a correlation matrix with no clusters,
#' "non-overlapping" for a correlation matrix with one layer of clusters, and
#' "overlapping" for a correlation matrix with up to 4 layers and a set number
#' of clusters per layer.
#' @param num_clusters if clusters = "non-overlapping" or clusters = "none" then
#' num_clusters is an integer value specifying the number of clusters. If clusters =
#' "overlapping" then num_clusters must be a vector, arranged in descending order, of
#' length equal to num_layers specifying the number of clusters per layer.
#' @param num_layers an positive integer value between 1 and 4, specifying the number
#' of cluster layers. Only needed if using "overlapping" clusters.
#' @return this function returns a D by D correlation matrix.
#'
#' @import propagate
#'
#' @examples
#' \dontrun{
#' library(ggcorrplot)
#'
#' ### This generates a 50 by 50 correlation matrix with no clusters.
#' gen_corr(D = 50, clusters = "none) %>%
#'          ggcorrplot(title = "Overlapping clusters")
#'
#' ### This generates a 50 by 50 correlation matrix with 5 non-overlapping clusters.
#' gen_corr(D = 50, clusters = "non-overlapping) %>%
#'          ggcorrplot(title = "Overlapping clusters")
#'
#' ### This generates a 60 by 60 correlation matrix consisting
#' ### of 4 layers with 10, 5, 3 and 2 clusters respectively.
#' gen_corr(D = 60,
#'          clusters = "overlapping",
#'          num_layers = 4,
#'          num_clusters = c(10,5,3,2)) %>%
#'                   ggcorrplot(title = "Overlapping clusters")
#'
#' }
#' @export
gen_corr <- function (D = 50,
                      clusters = c("none", "non-overlapping", "overlapping"),
                      num_clusters = NULL,
                      num_layers = NULL) {

        Grps <- num_clusters
        #set.seed(123)

        if (!(clusters %in%  c("none", "non-overlapping", "overlapping"))) stop("Please provide a valid clusters argument")

        if(clusters == "none"){
                # Unclustered covariance matrix
                Sigma <- diag(D)
                for (i in 1:D) for (j in 1:D) Sigma[i,j] <- 0.9^abs(i-j)
                Sigma <- propagate::cor2cov(Sigma, runif(D, 1, 5))
                corr <- cov2cor(Sigma)
        } else

                if(clusters == "non-overlapping"){
                        #----------------------
                        # distinct non-overlapping clusters:
                        #----------------------

                        if(is.null(num_clusters)) stop("Please provide a valid num_clusters argument when using Overlapping clusters")


                        Sigma <- matrix(0.9, D, D)
                        diag(Sigma) <- 1


                        for (i in 1:Grps) {
                                ix <- seq((i-1) * D / Grps + 1, i * D / Grps)
                                Sigma[ix, -ix] <- 0.0001                       #think about
                        }
                        Sigma <- propagate::cor2cov(Sigma, runif(D, 1, 5))
                        corr <- cov2cor(Sigma)
                } else

                        if(clusters == "overlapping"){
                                #----------------------
                                # distinct overlapping clusters:
                                #----------------------

                                if(is.null(num_layers)|num_layers<2){
                                        stop("Please provide a valid num_layers argument when using Overlapping clusters")
                                }else
                                        if(length(num_clusters) != num_layers){
                                                stop("Please provide a num_clusters argument with length equal to num_layers")
                                        }


                                Sigma <- matrix(0.6, D, D)
                                diag(Sigma) <- 1

                                for (i in 1:Grps[1]) {
                                        ix <- seq((i-1) * D / Grps[1] + 1, i * D / Grps[1])
                                        Sigma[ix, -ix] <- 0.7
                                }
                                if(num_layers>=2){
                                        for (i in 1:Grps[2]) {
                                                ix <- seq((i-1) * D / Grps[2] + 1, i * D / Grps[2])
                                                Sigma[ix, -ix] <- 0.5
                                        } }
                                if(num_layers>=3){
                                        for (i in 1:Grps[3]) {
                                                ix <- seq((i-1) * D / Grps[3] + 1, i * D / Grps[3])
                                                Sigma[ix, -ix] <- 0.3
                                        } }
                                if(num_layers>=4){
                                        for (i in 1:Grps[4]) {
                                                ix <- seq((i-1) * D / Grps[4] + 1, i * D / Grps[4])
                                                Sigma[ix, -ix] <- 0.05
                                        } }
                        }

        Sigma <- propagate::cor2cov(Sigma, runif(D, 1, 5))  #Is this necessary???
        corr <- cov2cor(Sigma)

        return(corr)

}
