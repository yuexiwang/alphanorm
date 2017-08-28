#####
# Introduction of the package
#####

#' alphanorm: A package for alpha-norm regularization model
#'
#' @description This package fits the alpha-norm regularization path for regression via
#' cyclic coordinate descent and o proximal operator. It is useful in extra sparse and
#' highly correlated model.
#'
#' @details  The alphanorm package provides five function:
#' \code{alphanorm}, \code{coef.alphanorm}, \code{cv.alphanorm}, \code{plot.alphanorm}
#'  and \code{predict.alphanorm}
#'
#'It accepts x and y for regression model and is very flexible in the choice of tuning
#'pararmeters q and lambda. \code{cv.alphanorm} can help select the best tuning parameters
#'using cross-validation. \code{plot.alphanorm} can produce the regularization path over
#'a grid of values for lambda.
#'
#'@author Guanhao Feng, Nicholas G Polson, Yuexi Wang and Jianeng Xu
#'
#'Maintainer: Yuexi Wang <yxwang99@uchicago.edu>
#'
#' @references
#' Feng, Guanhao and Polson, Nicholas G and Wang, Yuexi and Xu, Jianeng,
#' Sparse Regularization in Marketing and Economics (August 20, 2017).
#' Available at SSRN: \url{https://ssrn.com/abstract=3022856}
#'
#' Marjanovic, G. and V. Solo (2014). lq sparsity penalized linear regression with cyclic descent.
#' IEEE Transactions on Signal Processing 62(6), 1464–1475.
#'
#'
#' @docType package
#' @name alphanorm-package
NULL
