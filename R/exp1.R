.reshapeexp1 <- function(x, scale, matrix = FALSE) {
    
   n <- length(x)
    
    if (n > 1L) {
        msg <- paste("when 'x' has length n > 1, 'scale' and 'shape'",
                     " must have length 'n' or 1")
        if (length(scale) != n) {
            if (length(scale) != 1L) stop(msg)
            scale <- rep(scale, n)
        }
        nms <- names(x)
    } else {
        
        n <- length(scale)
        
        if (n > 1L) {
            x <- rep(x, n)
            scale <- rep(scale, length.out = n)
        }
        
    }
    
   if (matrix)  {
       
       return(cbind("x" = x, "scale" = scale))
       
   } else {
       
       if (length(x) > 1L) nms <- names(x)
       else if (length(scale) > 1L) nms <- names(scale)
       else nms <- ""
       
       return(list("n" = n, "x" = x,
                   "scale" = scale, 
                   "nms" = nms))
   }
}

## ***********************************************************************
##
##' @description Density, distribution function, quantile function and random
##' generation for the one-parameter Exponential Distribution
##' distribution with scale parameter \code{scale}.
##'
##' @details The survival and density functions are given by
##' \deqn{S(x) = \exp\{-x / \sigma\} \qquad
##'       f(x) = \frac{1}{\sigma} \exp\{-x / \sigma\} \qquad (x > 0)}{
##'       S(x) = exp(-x / sigma)  f(x) = exp(-x / sigma) / sigma (x > 0)}
##' where \eqn{\sigma} is the scale parameter. This distribution is
##' the Generalized Pareto Distribution for a shape \eqn{\xi = 0}.
##'
##' @name Exp1
##' @rdname Exp1
##' @export
##' 
##' @aliases pexp1 qexp1 rexp1
##' 
##' @title Density, Distribution Function, Quantile Function and
##'     Random Generation for the One-Parameter Exponential
##'     Distribution
##' 
##' @param scale Scale parameter. Numeric vector with suitable length,
##'     see \bold{Details}. Can not contain non-finite value.
##'
##' @param log Logical; if \code{TRUE}, densities \code{p} are
##'     returned as \code{log(p)}.
##' 
##' @param deriv Logical. If \code{TRUE}, the gradient of each
##'     computed value w.r.t. the parameter vector is computed, and
##'     returned as a \code{"gradient"} attribute of the result. This
##'     is a numeric array with dimension \code{c(n, 1)} where
##'     \code{n} is the length of the first argument, i.e. \code{x},
##'     \code{p} or \code{q}, depending on the function.
##'
##' @param hessian Logical. If \code{TRUE}, the Hessian of each
##'     computed value w.r.t. the parameter vector is computed, and
##'     returned as a \code{"hessian"} attribute of the result. This
##'     is a numeric array with dimension \code{c(n, 1, 1)} where
##'     \code{n} is the length of the first argument, i.e. \code{x},
##'     \code{p} or depending on the function.
##' 
##' @param array Logical. If \code{TRUE}, the simulated values form a
##'     numeric matrix with \code{n} columns and \code{np} rows where
##'     \code{np} is the number of exponential parameter values i.e.,
##'     the length of \code{scale}. This option is useful to cope with
##'     so-called \emph{non-stationary} models with exponential
##'     margins. See \bold{Examples}. The default value is
##'     \code{length(scale) > 1}.
##'
##' @param x,q Vector of quantiles.
##'
##' @param p Vector of probabilities.
##'
##' @param n Sample size.
##' 
##' @param lower.tail Logical; if \code{TRUE} (default), probabilities
##'     are \eqn{\textrm{Pr}[X \leq x]}{Pr[X <= x]}, otherwise,
##'     \eqn{\textrm{Pr}[X > x]}{Pr[X < x]}.
##'
##' @return A numeric vector with its length equal to the maximum of
##'     the two lengths: that of the first argument and that of the
##'     parameter \code{scale}. When \code{deriv} is \code{TRUE}, the
##'     returned value has an attribute named \code{"gradient"} which
##'     is a matrix with \eqn{n} lines and \eqn{1} column containing
##'     the derivative. A row contains the partial derivative of the
##'     corresponding element w.r.t. the parameter \code{"scale"}.
##'
##' @details The probability functions \code{d}, \code{p} and \code{q}
##'     all allow the parameter \code{scale} to be a vector. Then the
##'     recycling rule is used to get two vectors of the same length,
##'     corresponding to the first argument and to the scale
##'     parameter. This behaviour is the standard one for the
##'     probability functions of the \strong{stats} package but is
##'     unusual in R packages devoted to Extreme Value in which the
##'     parameters must generally have length one. Note that the
##'     provided functions can be used e.g. to evaluate the quantile
##'     with a given probability for a large number of values of the
##'     parameter vector \code{shape}. This is frequently required in
##'     he Bayesian framework with MCMC inference.
##'
##' @note The attributes \code{"gradient"} and \code{"hessian"} have
##'     dimension \code{c(n, 1)} and \code{c(n, 1, 1)}, even when
##'     \code{n} equals \code{1}. Use the \code{drop} method on these
##'     objects to drop the extra dimension if wanted i.e. to get a
##'     gradient vector and a Hessian matrix.
##'
##' @seealso The exponential distribution
##'     \code{\link[stats]{Exponential}} with \eqn{rate} being the
##'     inverse scale.
##' 
##' @examples
##' ## Illustrate the effect of recycling rule.
##' pexp1(1.0, scale = 1:4, lower.tail = FALSE) - exp(-1.0 / (1:4))
##' pexp1(1:4, scale = 1:4, lower.tail = FALSE) - exp(-1.0)
##'
##' ## With gradient and Hessian.
##' pexp1(c(1.1, 1.7), scale = 1, deriv = TRUE, hessian = TRUE)
##'
##' ti <- 1:60; names(ti) <- 2000 + ti
##' sigma <- 1.0 + 0.7 * ti
##' ## simulate 40 paths
##' y <- rexp1(n = 40, scale = sigma)
##' matplot(ti, y, type = "l", col = "gray", main = "varying scale")
##' lines(ti, apply(y, 1, mean))
##' 
dexp1 <- function(x, scale = 1.0, log = FALSE,
                  deriv = FALSE, hessian = FALSE) {
    
    if (!all(is.finite(scale))) {
        stop("exp1 parameter must be finite (non NA)")
    }
    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'gradient' is TRUE")
    }
    
    res <- .Call(Call_dexp1,
                 as.double(x),
                 as.double(scale),
                 as.integer(log),
                 as.integer(deriv),
                 as.integer(hessian))

    names(res) <- names(x)
  
    if (deriv) {
        n <- length(res)
        nm1 <- c("scale")
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 1L),
                  dimnames = list(names(x), nm1))
        
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 1L, 1L),
                      dimnames = list(names(x), nm1, nm1))
        }
        
    }

    return(res)
    
}

##' @name Exp1
##' @rdname Exp1
##' @export
pexp1 <- function(q, scale = 1.0, lower.tail = TRUE,
                  deriv = FALSE, hessian = FALSE) {

    if (!all(is.finite(scale))) {
        stop("exp1 parameter must be finite (non NA)")
    }
    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'gradient' is TRUE")
    }

    res <- .Call(Call_pexp1,
                 as.double(q),
                 as.double(scale),
                 as.integer(lower.tail),
                 as.integer(deriv),
                 as.integer(hessian))
    
    names(res) <- names(q)
    
    if (deriv) {
        n <- length(res)
        nm1 <- c("scale")
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 1L),
                  dimnames = list(names(q), c("scale")))
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 1L, 1L),
                      dimnames = list(names(q), nm1, nm1))
        }
    }
    return(res)
    
}

##' @name Exp1
##' @rdname Exp1
##' @export
qexp1 <- function(p, scale = 1.0, lower.tail = TRUE,
                  deriv = FALSE, hessian = FALSE) {

    if (!all(is.finite(scale))) {
        stop("exp1 parameter must be finite (non NA)")
    }
    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'gradient' is TRUE")
    }
    
    if (min(p, na.rm = TRUE) < 0.0 || max(p, na.rm = TRUE) > 1.0) 
        stop("`p' must contain probabilities in [0, 1]")
    
    res <- .Call(Call_qexp1,
                 as.double(p),
                 as.double(scale),
                 as.integer(lower.tail),
                 as.integer(deriv),
                 as.integer(hessian))
    
    names(res) <- names(p)
    
    if (deriv) {
        n <- length(res)
        nm1 <- c("scale")
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 1L),
                  dimnames = list(names(p), nm1))
        
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 1L, 1L),
                      dimnames = list(names(p), nm1, nm1))
        }
    }
    return(res)   

}

##' @name Exp1
##' @rdname Exp1
##' @export
rexp1 <- function(n, scale = 1.0, array) {
    
    if (!all(is.finite(scale))) {
        stop("exp1 parameter must be finite (non NA)")
    }

    if (missing(array)) {
        array <- length(scale) > 1L
    }
    
    if (array) {
        
        ## If 'array' is TRUE we return a matrix having the simulated values
        ## as its columns
        L <- .reshapeexp1(x = 1.0, scale = scale)
        
        res <- array(NA, dim = c(L$n, n),
                     dimnames = list(L$nms, paste("sim", 1:n, sep = "")) )
        
        scale <- array(L$scale, dim = c(L$n, n))
        
        nn <- L$n * n
        res <- array(qexp1(runif(nn), scale = scale),
                     dim = c(L$n, n),
                     dimnames = list(NULL, paste("sim", 1:n, sep = "")))
        
    } else {
        res <- qexp1(runif(n), scale = scale)
    }

    res 
}

