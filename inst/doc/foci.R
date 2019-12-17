## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, results='hide', message=FALSE, warning=FALSE----------------------
library(FOCI)

## -----------------------------------------------------------------------------
n = 1000
p = 150
X = matrix(rnorm(n * p), ncol = p)
Y = X[, 1] * X[, 2] + sin(X[, 1] * X[, 3]) + X[, 4]^2
result = foci(Y, X)
result

