##'
##' @title \packageTitle{nieve}
##' 
##' @description
##' 
##' The DESCRIPTION file:
##' \packageDESCRIPTION{nieve}
##' \packageIndices{nieve}
##' 
##' The \pkg{nieve} package provides utility functions for Extreme
##' Value Analysis. It includes the probability functions for the
##' two-parameter Generalized Pareto Distribution (GPD) and for the
##' three-parameter Generalized Extreme Value (GEV)
##' distribution. These functions are vectorized w.r.t. the parameters
##' and optionnaly provide the exact derivatives w.r.t. the
##' parameters: gradient and Hessian which can be used in optimisation
##' e.g., to maximise the log-likelihood. Since the gradient is
##' available for the distribution function, the exact gradient of the
##' log-likelihood function is available even when censored
##' observations are used.
##'
##' These functions should behave like the probability functions of
##' the \pkg{stats} package: when a probability \code{p = 0.0} or
##' \code{p = 1.0} is given, the quantile functions should return the
##' lower and the upper end-point, be they finite or not. Also when
##' evaluated at \code{-Inf} and \code{Inf} the probability functions
##' should return \code{0.0} and \code{1.0}.
##' 
##' @docType package
##' @name nieve-package
##' @useDynLib nieve, .registration=TRUE
NULL
#> NULL
