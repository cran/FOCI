## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, results='hide', message=FALSE, warning=FALSE----------------------
library(FOCI)

## -----------------------------------------------------------------------------
n = 2000
p = 100
X = matrix(rnorm(n * p), ncol = p)
colnames(X) = paste0(rep("X", p), seq(1, p))
Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3]) + X[, 4]^2
result1 = foci(Y, X, numCores = 1)
result1

## -----------------------------------------------------------------------------
result2 = foci(Y, X, num_features = 5, stop = FALSE, numCores = 1)

