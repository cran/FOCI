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
x1 = matrix(runif(n), ncol = 1)
x2 = matrix(runif(n), ncol = 1)
x3 = matrix(runif(n), ncol = 1)
y = (x1 + x2 + x3) %% 1
# y is independent of each of x1 and x2 and x3 
codec(y, x1)
codec(y, x2)
codec(y, x3)

# y is independent of the pair (x1, x2)
codec(y, cbind(x1, x2))

# y is a function of (x1, x2, x3)
codec(y, cbind(x1, x2, x3))

# conditional on x3, y is a function of (x1, x2)
codec(y, cbind(x1, x2), x3)
# conditional on x3, y is independent of x1
codec(y, x1, x3)


