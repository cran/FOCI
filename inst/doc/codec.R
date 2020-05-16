## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_chunk$set(fig.width=7, fig.height=5) 

## ----setup, results='hide', message=FALSE, warning=FALSE----------------------
library(FOCI)

## -----------------------------------------------------------------------------
n = 10000
p = 3
x = matrix(runif(n * p), ncol = p)
y = (x[, 1] + x[, 2] + x[, 3]) %% 1
# y is independent of each of column of x 
codec(y, x[, 1])
codec(y, x[, 2])
codec(y, x[, 3])

# y is independent of the first two columns of x, x[, c(1, 2)]
codec(y, x[, c(1, 2)])

# y is a function of x
codec(y, x)

# conditional on the last column of x, y is a function of the first two columns
codec(y, x[, c(1, 2)], x[, 3])
# conditional on x[, 3], y is independent of x[, 1]
codec(y, x[, 1], x[, 3])


## -----------------------------------------------------------------------------
n = 1000
p = 2
x = matrix(rnorm(n * p), ncol = p)
y = x[, 1]^2 + x[, 2]^2
z = atan(x[, 1] / x[, 2])
# y is independent of z
codec(y, z)
# conditional on x[, 1], y is a function of z
codec(y, z, x[, 1])

