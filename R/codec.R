# codec -------------------------------------------------------------------------
#' Estimate the conditional dependence coefficient (CODEC)
#'
#' The conditional dependence coefficient (CODEC) is a measure of the amount of conditional dependence between
#' a random variable Y and a random vector Z given a random vector X, based on an i.i.d. sample of (Y, Z, X).
#' The coefficient is asymptotically guaranteed to be between 0 and 1.
#'
#' @param Y Vector (length n)
#' @param Z Matrix (n by q)
#' @param X Matrix (n by p), default is NULL
#' @param na.rm Remove NAs if TRUE
#' @details The value returned by codec can be positive or negative. Asymptotically, it is guaranteed
#' to be between 0 and 1. A small value indicates low conditional dependence between Y and Z given X, and
#' a high value indicates strong conditional dependence. The codec function is used by the  \code{\link{foci}} function
#' for variable selection.
#' @return The conditional dependence coefficient (CODEC) of Y and Z given X. If X == NULL, this is just a
#' measure of the dependence between Y and Z.
#' @import data.table
#' @importFrom stats complete.cases sd
#' @export
#' @author Mona Azadkia, Sourav Chatterjee, Norman Matloff
#' @references Azadkia, M. and Chatterjee, S. (2019). A simple measure
#' of conditional dependence.
#' \url{https://arxiv.org/pdf/1910.12327.pdf}.
#' @seealso \code{\link{foci}}, \code{\link[XICOR]{xicor}}
#' @examples
#' n = 1000
#' x <- matrix(runif(n * 2), nrow = n)
#' y <- (x[, 1] + x[, 2]) %% 1
#' # given x[, 1], y is a function of x[, 2]
#' codec(y, x[, 2], x[, 1])
#' # y is a function of x
#' codec(y, x)
#' z <- rnorm(n)
#' # y is a function of x given z
#' codec(y, x, z)
#' # y is independent of z given x
#' codec(y, z, x)
codec <- function(Y, Z, X = NULL, na.rm = TRUE){
  if(is.null(X)) {
    # if inputs are not in proper matrix format change if possible
    # otherwise send error
    if(!is.vector(Y)) {
      Y = as.vector(Y)
    }
    if(!is.matrix(Z)) {
      Z = as.matrix(Z)
    }
    if((length(Y) != nrow(Z))) stop("Number of rows of Y and X should be equal.")
    if (na.rm == TRUE) {
      # NAs are removed here:
      ok = complete.cases(Y,Z)
      Z = as.matrix(Z[ok,])
      Y = Y[ok]
    }

    n = length(Y)
    if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")

    q = ncol(Z)
    return(.estimateT(Y, Z))
  }
  # if inputs are not in proper matrix format change if possible
  # otherwise send error
  if(!is.vector(Y)) {
    Y = as.vector(Y)
  }
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  if(!is.matrix(Z)) {
    Z = as.matrix(Z)
  }
  if((length(Y) != nrow(X))) stop("Number of rows of Y and X should be equal.")
  if((length(Y) != nrow(Z))) stop("Number of rows of Y and Z should be equal.")
  if((nrow(Z) != nrow(X))) stop("Number of rows of Z and X should be equal.")
  if (na.rm == TRUE) {
    # NAs are removed here:
    ok = complete.cases(Y,Z,X)
    Z = as.matrix(Z[ok,])
    Y = Y[ok]
    X = as.matrix(X[ok,])
  }

  n = length(Y)
  if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")

  q = ncol(Z)
  p = ncol(X)

  return(.estimateConditionalT(Y, Z, X))
}
