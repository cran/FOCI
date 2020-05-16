#############################################
# HELPER FUNCTIONS for codec
#
# .estimateConditionalQ
# .estimateConditionalS
# .estimateConditionalT
# .estimateQ
# .estimateS
# .estimateT
##############################################


# .estimateConditionalQ -------------------------------------------------------------------------
# Estimate Q(Y, Z | X)
#
# Estimate Q(Y, Z | X), the numinator of the measure of conditional dependence of Y on Z given X
#
# @param X: Matrix of predictors (n by p)
# @param Z: Matrix of predictors (n by q)
# @param Y: Vector (length n)
#
# @return estimation \eqn{Q(Y, Z|X)}
.estimateConditionalQ <- function (Y, X, Z) {

  id <- group <- rnn <- NULL

  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  if(!is.matrix(Z)) {
    Z = as.matrix(Z)
  }

  n = length(Y)

  W = cbind(X, Z)

  # compute the nearest neighbor of X
  nn_X = RANN::nn2(X, query = X, k = 3)
  nn_index_X = nn_X$nn.idx[, 2]
  # handling repeated data
  repeat_data = which(nn_X$nn.dists[, 2] == 0)

  df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
  df_X[, rnn := .randomNN(id), by = "group"]

  nn_index_X[repeat_data] = df_X$rnn
  # nearest neighbors with ties
  ties = which(nn_X$nn.dists[, 2] == nn_X$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)

  if(length(ties) > 0) {
    helper_ties <- function(a) {
      distances <- proxy::dist(matrix(X[a, ], ncol = ncol(X)), matrix(X[-a, ], ncol = ncol(X)))
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }

    nn_index_X[ties] = sapply(ties, helper_ties)
  }

  # compute the nearest neighbor of W
  nn_W = RANN::nn2(W, query = W, k = 3)
  nn_index_W = nn_W$nn.idx[, 2]
  repeat_data = which(nn_W$nn.dists[, 2] == 0)

  df_W = data.table::data.table(id = repeat_data, group = nn_W$nn.idx[repeat_data])
  df_W[, rnn := .randomNN(id), by = "group"]

  nn_index_W[repeat_data] = df_W$rnn
  # nearest neighbors with ties
  ties = which(nn_W$nn.dists[, 2] == nn_W$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)

  if(length(ties) > 0) {
    helper_ties <- function(a) {
      distances <- proxy::dist(matrix(X[a, ], ncol = ncol(X)), matrix(X[-a, ], ncol = ncol(X)))
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }

    nn_index_W[ties] = sapply(ties, helper_ties)
  }

  # estimate Q
  R_Y = rank(Y, ties.method = "max")
  Q_n = sum(apply(rbind(R_Y, R_Y[nn_index_W]), 2, min),
            -apply(rbind(R_Y, R_Y[nn_index_X]), 2, min)) / (n^2)
  return(Q_n)
}




# .estimateConditionalS -------------------------------------------------------------------------
# Estimate S(Y, X)
#
# Estimate S(Y, X), the denuminator of the measure of dependence of Y on Z given X
#
# @param X: Matrix of predictors (n by p)
# @param Y: Vector (length n)
#
# @return estimation \eqn{S(Y, X)}
.estimateConditionalS <- function (Y, X){

  id <- group <- rnn <- NULL

  if(!is.matrix(X)) {
    X = as.matrix(X)
  }
  n = length(Y)

  # compute the nearest neighbor of X
  nn_X = RANN::nn2(X, query = X, k = 3)
  nn_index_X = nn_X$nn.idx[, 2]
  repeat_data = which(nn_X$nn.dists[, 2] == 0)

  df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
  df_X[, rnn := .randomNN(id), by = "group"]

  nn_index_X[repeat_data] = df_X$rnn
  # nearest neighbors with ties
  ties = which(nn_X$nn.dists[, 2] == nn_X$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)

  if(length(ties) > 0) {
    helper_ties <- function(a) {
      distances <- proxy::dist(matrix(X[a, ], ncol = ncol(X)), matrix(X[-a, ], ncol = ncol(X)))
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }

    nn_index_X[ties] = sapply(ties, helper_ties)
  }

  # estimate S
  R_Y = rank(Y, ties.method = "max")
  S_n = sum(R_Y - apply(rbind(R_Y, R_Y[nn_index_X]), 2, min)) / (n^2)

  return(S_n)
}


# estimateConditionalT -------------------------------------------------------------------------
# Estimate T(Y, Z | X)
#
# Estimate T(Y, Z | X), the measure of dependence of Y on Z given X
#
# @param Y: Vector (length n)
# @param Z: Matrix of predictors (n by q)
# @param X: Matrix of predictors (n by p)
#
# @return estimation of \eqn{T(Y, Z|X)}.
.estimateConditionalT <- function(Y, Z, X){
  S = .estimateConditionalS(Y, X)

  # happens only if Y is constant
  if (S == 0) {
    return(1)
  } else {
    return(.estimateConditionalQ(Y, X, Z) / S)
  }
}




# .estimateQ -------------------------------------------------------------------------
# Estimate Q(Y, X)
#
# Estimate Q(Y, X), the numinator of the measure of dependence of Y on X
#
# @param X: Matrix of predictors (n by p).
# @param Y: Vector (length n).
#
# @return estimation of \eqn{Q(Y, X)}.
.estimateQ <- function(Y, X) {

  id <- group <- rnn <- NULL

  if(!is.matrix(X)) {
    X = as.matrix(X)
  }

  n = length(Y)
  nn_X = RANN::nn2(X, query = X, k = 3)
  # remove the first nearest neighbor for each x which is x itself in case of no repeat data
  # when there is repeated data this is wrong but for repeated data we find the nearest
  # neighbors separately.
  nn_index_X = nn_X$nn.idx[, 2]

  # find all data points that are not unique
  repeat_data = which(nn_X$nn.dists[, 2] == 0)

  # for the repeated data points, choose one of their identicals at random and set its index
  # as the index of the nearest neighbor
  df_X = data.table::data.table(id = repeat_data, group = nn_X$nn.idx[repeat_data, 1])
  df_X[, rnn := .randomNN(id), by = "group"]
  nn_index_X[repeat_data] = df_X$rnn

  # nearest neighbors with ties
  ties = which(nn_X$nn.dists[, 2] == nn_X$nn.dists[, 3])
  ties = setdiff(ties, repeat_data)
  if(length(ties) > 0) {
    helper_ties <- function(a) {
      distances <- proxy::dist(matrix(X[a, ], ncol = ncol(X)), matrix(X[-a, ], ncol = ncol(X)))
      ids <- which(distances == min(distances))
      x <- sample(ids, 1)
      return(x + (x >= a))
    }

    nn_index_X[ties] = sapply(ties, helper_ties)
  }


  # estimate Q
  R_Y = rank(Y, ties.method = "max")
  L_Y = rank(-Y, ties.method = "max")
  Q_n = sum(apply(rbind(R_Y, R_Y[nn_index_X]), 2, min) - (L_Y^2)/n) / (n^2)

  return(Q_n)
}



# .estimateS -------------------------------------------------------------------------
# Estimate S(Y)
#
# Estimate S(Y) , the denuminator of the measure of dependence of Y on X
#
# @param Y: Vector (length n).
# @return estimation of \eqn{S(Y)}.
.estimateS <- function (Y) {
  n = length(Y)
  L_Y = rank(-Y, ties.method = "max")
  S_n = gmp::asNumeric(sum(gmp::as.bigz(L_Y) * gmp::as.bigz(n - L_Y))) / (n^3)
  return(S_n)
}



# .estimateT -------------------------------------------------------------------------
# Estimate T(Y, X)
#
# Estimate T(Y, X), the measure of dependence of Y on X
#
# @param X: Matrix of predictors (n by p)
# @param Y: Vector (length n)
# @return estimation of \eqn{T(Y, X) = Q(Y, X) / S(Y)}.
.estimateT <- function(Y, X){
  # May 15, Mona removed FOCI:::
  S = .estimateS(Y)
  # happens only if Y is a constant vector.
  if (S == 0) {
    return(1)
  } else {
    return(.estimateQ(Y, X) / S)
  }
}
