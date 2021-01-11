#' Correlation Matrices
#' @description This is a data set containing 30 empirical correlation matrices
#' obtained by randomly selecting a set of 50 stocks from the S&P 500 over a random
#' time period of 252 days between 1 January 2000 and 26 December 2020. The correlation
#' matrices were then labelled as ‘stressed’ if the equi-weighted basket of stocks had
#' a Sharpe ratio below -0.5, ‘rally’ if the equi-weighted basket of stocks had a Sharpe
#' ratio above 2  and ‘normal’ if the equi-weighted basket of stocks had a Sharpe
#' in-between -0.5 and 2.
#' @name corr_mats
#' @docType data
#' @usage data(corr_mats)
#' @author Nathan Potgieter  \email{nathan.potieter56@gmail.com}
#' @keywords data
#' @format A list containing three types of correlation matrices. They are
#' labeled corr_normal, corr_stressed and corr_normal. Each type contains 10
#' correlation matrices.
"corr_mats"
