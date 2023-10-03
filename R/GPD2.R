.reshapeGPD2 <- function(x, scale, shape, matrix = FALSE) {
   
   n <- length(x)
   
   if (n > 1L) {
       msg <- paste("when 'x' has length n > 1, 'scale' and 'shape'",
                    " must have length 'n' or 1")
       if (length(scale) != n) {
           if (length(scale) != 1L) stop(msg)
           scale <- rep(scale, n)
       }
       if (length(shape) != n) {
           if (length(shape) != 1L) stop(msg)
           shape <- rep(shape, n)
       }
       nms <- names(x)
   } else {

       n <- max(c(length(scale), length(shape)))
       
       if (n > 1L) {
           
           if (((length(scale) > 1L) && (length(scale) < n)) ||
               ((length(shape) > 1L) && (length(shape) < n)) ) {
               stop("when 'x' has length 1, the lengths of the non-scalar ",
                    "elements among 'scale' and 'shape' must be ",
                    "the same")
          }
           x <- rep(x, n)
           scale <- rep(scale, length.out = n)
           shape <- rep(shape, length.out = n) 
      }
       
   }
   
   if (matrix)  {

       return(cbind("x" = x, "scale" = scale, "shape" = shape))
      
   } else {

       if (length(x) > 1L) nms <- names(x)
       else if (length(scale) > 1L) nms <- names(scale)
       else if (length(shape) > 1L) nms <- names(shape)
       else nms <- ""
       
       return(list("n" = n, "x" = x,
                   "scale" = scale, "shape" = shape,
                   "nms" = nms))
   }
}

## ***********************************************************************

##' @description Density, distribution function, quantile function and
##'     random generation for the two-parameter Generalized Pareto
##'     Distribution (GPD) distribution with \code{scale} and
##'     \code{shape}.
##'
##' @name GPD2
##' @rdname GPD2
##' @export
##' 
##' @aliases pGPD2 qGPD2 rGPD2
##' 
##' @title Density, Distribution Function, Quantile Function and
##'     Random Generation for the Two-Parameter Generalized Pareto
##'     Distribution (GPD)
##' 
##' @param scale Scale parameter. Numeric vector with suitable length,
##'     see \bold{Details}.
##'
##' @param shape Shape parameter. Numeric vector with suitable length,
##'     see \bold{Details}.
##'
##' @param log Logical; if \code{TRUE}, densities \code{p} are
##'     returned as \code{log(p)}.
##' 
##' @param deriv Logical. If \code{TRUE}, the gradient of each
##'     computed value w.r.t. the parameter vector is computed, and
##'     returned as a \code{"gradient"} attribute of the result. This
##'     is a numeric array with dimension \code{c(n, 2)} where
##'     \code{n} is the length of the first argument, i.e. \code{x},
##'     \code{p} or \code{q}, depending on the function.
##'
##' @param hessian Logical. If \code{TRUE}, the Hessian of each
##'     computed value w.r.t. the parameter vector is computed, and
##'     returned as a \code{"hessian"} attribute of the result. This
##'     is a numeric array with dimension \code{c(n, 2, 2)} where
##'     \code{n} is the length of the first argument, i.e. \code{x},
##'     \code{p} or depending on the function.
##' 
##' @param array Logical. If \code{TRUE}, the simulated values form a
##'     numeric matrix with \code{n} columns and \code{np} rows where
##'     \code{np} is the number of GPD parameter values. This number
##'     is obtained by recycling the two GPD parameters vectors to a
##'     common length, so \code{np} is the maximum of the lengths of
##'     the parameter vectors \code{scale} and \code{shape}. This
##'     option is useful to cope with so-called \emph{non-stationary}
##'     models with GPD margins. See \bold{Examples}. The default
##'     value is \code{TRUE} if any of the vectors \code{scale} and
##'     \code{shape} has length \code{> 1} and \code{FALSE} otherwise.
##'
##' @param x,q Vector of quantiles.
##'
##' @param p Vector of probabilities.
##'
##' @param n Sample size.
##' 
##' @param lower.tail Logical; if \code{TRUE} (default), probabilities
##'     are \eqn{\textrm{Pr}[X \leq x]}{Pr[X <= x]}, otherwise,
##'     \eqn{\textrm{Pr}[X > x]}{Pr[X > x]}.
##'
##' @return A numeric vector with length equal to the maximum of the
##' four lengths: that of the first argument and that of the two
##' parameters \code{scale} and \code{shape}. When \code{deriv} is
##' \code{TRUE}, the returned value has an attribute named
##' \code{"gradient"} which is a matrix with \eqn{n} lines and \eqn{2}
##' columns containing the derivatives. A row contains the partial
##' derivatives of the corresponding element w.r.t. the two parameters
##' \code{"scale"} and \code{"shape"} in that order.
##'
##' @details
##' Let \eqn{\sigma >0} and \eqn{\xi} denote the scale and the shape; the
##' survival function \eqn{S(x) := \textrm{Pr}[X > x]}{Pr[X < x]} is given
##' for \eqn{x \geq 0}{x >= 0} by
##' \deqn{S(x) = \left[1 + \xi x/ \sigma \right]_+^{-1/\xi}}{[1 + \xi * x / \sigma]_+^(-1/\xi)} for \eqn{\xi \neq 0}{\xi != 0} where \eqn{v_+ := \max\{v, \, 0\}}{v_+ = max(v, 0)} 
##' and by \deqn{S(x) = \exp\{-x/\sigma\}}{exp(- x / \sigma)}
##' for \eqn{\xi = 0}{\xi = 0}. For \eqn{x < 0} we have \eqn{S(x) = 1}:
##' the support of the distribution is \eqn{(0,\,\infty(}{(0, Inf(}.
##'
##' The probability functions \code{d}, \code{p} and \code{q} all
##' allow each of the two GP parameters to be a vector. Then the
##' recycling rule is used to get three vectors of the same length,
##' corresponding to the first argument and to the two GP
##' parameters. This behaviour is the standard one for the probability
##' functions of the \strong{stats}. Note that the provided functions
##' can be used e.g.  to evaluate the quantile with a given
##' probability for a large number of values of the parameter vector
##' \code{c(shape, scale)}. This is frequently required in he Bayesian
##' framework with MCMC inference.
##'
##' @note The attributes \code{"gradient"} and \code{"hessian"} have
##' dimension \code{c(n, 2)} and \code{c(n, 2, 2)}, even when \code{n}
##' equals \code{1}. Use the \code{drop} method on these objects to
##' drop the extra dimension if wanted i.e. to get a gradient vector
##' and a Hessian matrix.
##' 
##' @examples
##' ## Illustrate the effect of recycling rule.
##' pGPD2(1.0, scale = 1:4, shape = 0.0, lower.tail = FALSE) - exp(-1.0 / (1:4))
##' pGPD2(1:4, scale = 1:4, shape = 0.0, lower.tail = FALSE) - exp(-1.0)
##'
##' ## With gradient and Hessian.
##' pGPD2(c(1.1, 1.7), scale = 1, shape = 0, deriv = TRUE, hessian = TRUE)
##'
##' ## simulate 40 paths
##' ti <- 1:20
##' names(ti) <- 2000 + ti
##' y <- rGPD2(n = 40, scale = ti, shape = 0.05)
##' matplot(ti, y, type = "l", col = "gray", main = "varying scale")
##' lines(ti, apply(y, 1, mean))
##' 
dGPD2 <- function(x, scale = 1.0, shape = 0.0,
                  log = FALSE,
                  deriv = FALSE, hessian = FALSE) {

    ## if (!all(is.finite(scale)) || !all(is.finite(shape))) {
    ##     stop("GPD2 parameters must be finite (non NA)")
    ## }
        
    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'gradient' is TRUE")
    }
    
    res <- .Call(Call_dGPD2,
                 as.double(x),
                 as.double(scale),
                 as.double(shape),
                 as.integer(log),
                 as.integer(deriv),
                 as.integer(hessian))
    
    names(res) <- names(x)
    
    if (deriv) {
        n <- length(res)
        nm2 <- c("scale", "shape")
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 2L),
                  dimnames = list(names(x), nm2))
        
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 2L, 2L),
                      dimnames = list(names(x), nm2, nm2))
        }
        
    }

    return(res)
    
}

##' @rdname GPD2
##' @export
pGPD2 <- function(q, scale = 1.0, shape = 0.0, lower.tail = TRUE,
                  deriv = FALSE, hessian = FALSE) {

    ## if (!all(is.finite(scale)) || !all(is.finite(shape))) {
    ##     stop("GPD2 parameters must be finite (non NA)")
    ## }
    
    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'gradient' is TRUE")
    }

    res <- .Call(Call_pGPD2,
                 as.double(q),
                 as.double(scale),
                 as.double(shape),
                 as.integer(lower.tail),
                 as.integer(deriv),
                 as.integer(hessian))
    
    names(res) <- names(q)
  
    if (deriv) {
        n <- length(res)
        nm2 <- c("scale", "shape")
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 2L),
                  dimnames = list(names(q), c("scale", "shape")))
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 2L, 2L),
                      dimnames = list(names(q), nm2, nm2))
        }
    }
    return(res)
    
}

##' @rdname GPD2
##' @export
qGPD2 <- function(p, scale = 1.0, shape = 0.0, lower.tail = TRUE,
                  deriv = FALSE, hessian = FALSE) {
    
    ## if (!all(is.finite(scale)) || !all(is.finite(shape))) {
    ##     stop("GPD2 parameters must be finite (non NA)")
    ## }

    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'gradient' is TRUE")
    }
    
    if (min(p, na.rm = TRUE) < 0.0 || max(p, na.rm = TRUE) > 1.0) 
        stop("`p' must contain probabilities in [0, 1]")
    
    res <- .Call(Call_qGPD2,
                 as.double(p),
                 as.double(scale),
                 as.double(shape),
                 as.integer(lower.tail),
                 as.integer(deriv),
                 as.integer(hessian))

    names(res) <- names(p)
    
    if (deriv) {
        n <- length(res)
        nm2 <- c("scale", "shape")
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 2L),
                  dimnames = list(names(p), nm2))
        
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 2L, 2L),
                      dimnames = list(names(p), nm2, nm2))
        }
    }
    return(res)   

}

##' @rdname GPD2
##' @importFrom stats runif
##' @export
rGPD2 <- function(n, scale = 1.0, shape = 0.0, array) {

    ## if (!all(is.finite(scale)) || !all(is.finite(shape))) {
    ##     stop("GPD2 parameters must be finite (non NA)")
    ## }

    if (missing(array)) {
        array <- max(c(length(scale), length(shape))) > 1L
    }
    
    if (array) {
        
        ## If 'array' is TRUE we return a matrix having the simulated values
        ## as its columns
        L <- .reshapeGPD2(x = 1.0, scale = scale, shape = shape)
        
        res <- array(NA, dim = c(L$n, n),
                     dimnames = list(L$nms, paste("sim", 1:n, sep = "")) )
        
        scale <- array(L$scale, dim = c(L$n, n))
        shape <- array(L$shape, dim = c(L$n, n))
        
        nn <- L$n * n
        res <- array(qGPD2(runif(nn), scale = scale, shape = shape),
                     dim = c(L$n, n),
                     dimnames = list(NULL, paste("sim", 1:n, sep = "")))
        
    } else {
        res <- qGPD2(runif(n), scale = scale, shape = shape)
    }

    res 
}

