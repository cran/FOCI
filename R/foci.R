#############################################
# MAIN FUNCTIONS:
# foci:  Performs feature ordering by conditional independence
##############################################

# foci -------------------------------------------------------------------------
#' Variable selection by the FOCI algorithm
#'
#' foci is a variable selection algorithm based on the measure of conditional dependence \code{\link{codec}}.
#'
#' @param X Matrix of predictors (n by p)
#' @param Y Vector of responses (length n)
#' @param num_features Number of variables to be selected, cannot be larger than p. The default value is NULL and in that
#' case it will be set equal to p. If stop == TRUE (see below), then num_features is irrelevant.
#' @param stop Stops at the first instance of negative codec, if TRUE.
#' @param na.rm Removes NAs if TRUE.
#' @param factor Converts factor variables to integers if TRUE.
#' @param standardize Standardize covariates if TRUE.
#' @details foci is a forward stepwise algorithm that uses the conditional dependence coefficient (\code{\link{codec}})
#' defined in Azadkla and Chatterjee (2019) at each step, instead of the multiple correlation coefficient
#' as in ordinary forward stepwise. If stop == TRUE, the process is stopped at the first instance of
#' nonpositive codec, thereby selecting subset of variables. Otherwise, a set of covariates of size
#' num_features, ordered according to predictive power (as measured by codec) is produced.
#' @return A vector of selected covariates, in decreasing order of predictive power.
#' @export
#' @author Mona Azadkia and Sourav Chatterjee
#' @references Azadkia, M. and Chatterjee, S. (2019). A simple measure of conditional dependence. \url{https://arxiv.org/pdf/1910.12327.pdf}
#' @seealso \code{\link{codec}}
#' @examples
#' # Example 1
#' n = 1000
#' p = 100
#' x <- matrix(rnorm(n * p), nrow = n)
#' y <- x[, 1] * x[, 10] + x[, 20]^2
#' # with num_features equal to 3 and stop equal to FALSE, foci will give a list of
#' # three selected features
#' result1 = foci(y, x, num_features = 3, stop = FALSE)
#' result1
#' # Example 2
#' # foci for the same example, this time it will stop according to its stopping rule.
#' result2 = foci(y, x)
#' result2
#'
foci <- function(Y, X, num_features = NULL, stop = TRUE, na.rm = TRUE, factor = TRUE, standardize = TRUE){
  # if inputs are not in proper format change if possible
  # otherwise send error
  if(!is.vector(Y)) {
    Y = as.vector(Y)
  }
  if(!is.matrix(X)) {
    X = as.matrix(X)
  }

  if (is.null(num_features)) num_features = dim(X)[2]

  if((length(Y) != nrow(X))) stop("Number of rows of Y and X should be equal.")

  if (na.rm == TRUE) {
    # NAs are removed here:
    ok = complete.cases(Y,X)
    X = as.matrix(X[ok,])
    Y = Y[ok]
  }


  n = length(Y)

  if(n < 2) stop("Number of rows with no NAs should be bigger than 1.")

  p = ncol(X)

  if (num_features > p) stop("Number of features should not be larger than maximum number of original features.")
  if ((floor(num_features) != num_features) || (num_features <= 0)) stop("Number of features should be a positive integer.")

  # Convert factor columns into numeric
  if (factor == TRUE) {
    if (!is.numeric(Y)) Y = as.numeric(factor(Y))
    X1 = matrix(nrow = n, ncol = p)
    for (i in 1:p) {
      a = X[,i]
      if (!is.numeric(a)) X1[,i] = as.numeric(factor(a))
      if (is.numeric(a)) X1[,i] = X[,i]
    }
    X = X1
  }

  if (standardize == TRUE) {
    for (i in 1:p) X[,i] = (X[,i] - mean(X[,i]))/sd(X[,i])
  }

  Q = rep(0, num_features)
  index_select = rep(0, num_features)
  # select the first variable
  estimateQFixedY <- function(x){
    return(.estimateQ(Y, x))
  }

  if (is.null(dim(X))) {
    seq_Q = estimateQFixedY(X)
  } else {
    seq_Q = apply(X, 2, estimateQFixedY)
  }

  Q[1] = max(seq_Q)
  if (Q[1] <= 0 & stop == TRUE) return(NULL)
  index_max = min(which(seq_Q == Q[1]))
  index_select[1] = index_max
  count = 1

  # select rest of the variables
  while (count < num_features) {
    seq_Q = rep(0, p - count)
    # indices that have not been selected yet
    index_left = setdiff(seq(1, p), index_select[1:count])

    # find the next best feature
    estimateQFixedYandSubX <- function(x){
      return(.estimateQ(Y, cbind(X[, index_select[1:count]], x)))
    }

    if (length(index_left) == 1) {
      seq_Q = estimateQFixedYandSubX(X[, index_left])
    } else {
      seq_Q = apply(X[, index_left], 2, estimateQFixedYandSubX)
    }
    Q[count + 1] = max(seq_Q)
    index_max = min(which(seq_Q == Q[count + 1]))
    if (Q[count + 1] <= Q[count] & stop == TRUE) break
    index_select[count + 1] = index_left[index_max]
    count = count + 1
  }

  return(index_select[1:count])
}
