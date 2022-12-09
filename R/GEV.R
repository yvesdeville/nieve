.reshapeGEV <- function(x, loc, scale, shape, matrix = FALSE) {
   
   n <- length(x)
   
   if (n > 1L) {
      msg <- paste("when 'x' has length n > 1, 'loc', 'scale' and 'shape'",
                   " must have length 'n' or 1")
      if (length(loc) != n) {
         if (length(loc) != 1L) stop(msg)
         loc <- rep(loc, n)
      }
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

       n <- max(c(length(loc), length(scale), length(shape)))
       
       if (n > 1L) {
          
          if (((length(loc) > 1L) && (length(loc) < n)) ||
              ((length(scale) > 1L) && (length(scale) < n)) ||
              ((length(shape) > 1L) && (length(shape) < n)) ) {
              stop("when 'x' has length 1, the lengths of the non-scalar ",
                   "elements among 'loc', 'scale' and 'shape' must be ",
                   "the same")
          }
          x <- rep(x, n)
          loc <- rep(loc, length.out = n)
          scale <- rep(scale, length.out = n)
          shape <- rep(shape, length.out = n) 
      }
       
   }
   
   if (matrix)  {

       return(cbind("x" = x, "loc" = loc, "scale" = scale, "shape" = shape))
      
   } else {

       if (length(x) > 1L) nms <- names(x)
       else if (length(loc) > 1L) nms <- names(loc)
       else if (length(scale) > 1L) nms <- names(scale)
       else if (length(shape) > 1L) nms <- names(shape)
       else nms <- ""
       
       return(list("n" = n, "x" = x,
                   "loc" = loc, "scale" = scale, "shape" = shape,
                   "nms" = nms))
   }
}

##  ****************************************************************************
##' @description Density, distribution function, quantile function and
##'     random generation for the Generalized Extreme Value (GEV)
##'     distribution with parameters \code{loc}, \code{scale} and
##'     \code{shape}.
##'     The distribution function \eqn{F(x) = \textrm{Pr}[X \leq x]}{F(x)= Pr[X <= x]}
##'     is given by 
##'     \deqn{F(x) = \exp\left\{-[1 + \xi z]^{-1/\xi}\right\}}{F(x) = exp(-(1 + xi * z))}
##'     when \eqn{\xi \neq 0}{xi != 0} and \eqn{1 + \xi z > 0}{1 + xi * z > 0}, and by 
##'     \deqn{F(x) = \exp\left\{-e^{-z}\right\}}{F(x) = exp(-exp(-z))}
##'     for \eqn{\xi =0}{xi == 0} where \eqn{z := (x - \mu) / \sigma} in both cases.
##' @name GEV
##' @rdname GEV
##' @export
##' 
##' @title Density, Distribution Function, Quantile Function and
##' Random Generation for the Generalized Extreme Value (GEV)
##' Distribution
##'
##' @param loc Location parameter. Numeric vector with suitable
##' length, see \bold{Details}. Can not contain no-finite value.
##'
##' @param scale Scale parameter. Numeric vector with suitable length,
##' see \bold{Details}. Can not contain no-finite value.
##'
##' @param shape Shape parameter. Numeric vector with suitable length,
##' see \bold{Details}. Can not contain non-finite value.
##'
##' @param log Logical; if \code{TRUE}, densities \code{p} are
##' returned as \code{log(p)}.
##' 
##' @param deriv Logical. If \code{TRUE}, the gradient of each
##'     computed value w.r.t. the parameter vector is computed, and
##'     returned as a \code{"gradient"} attribute of the result. This
##'     is a numeric array with dimension \code{c(n, 3)} where
##'     \code{n} is the length of the first argument, i.e. \code{x},
##'     \code{p} or \code{q} depending on the function.
##'
##' @param hessian Logical. If \code{TRUE}, the Hessian of each
##'     computed value w.r.t. the parameter vector is computed, and
##'     returned as a \code{"hessian"} attribute of the result. This
##'     is a numeric array with dimension \code{c(n, 3, 3)} where
##'     \code{n} is the length of the first argument, i.e. \code{x},
##'     \code{p} or depending on the function.
##'
##' @param array Logical. If \code{TRUE}, the simulated values form a
##'     numeric matrix with \code{n} columns and \code{np} rows where
##'     \code{np} is the number of GEV parameter values. This number
##'     is obtained by recycling the three GEV parameters vectors to a
##'     common length, so \code{np} is the maximum of the lengths of
##'     the parameter vectors \code{loc}, \code{scale},
##'     \code{shape}. This option is useful to cope with so-called
##'     \emph{non-stationary} models with GEV margins. See
##'     \bold{Examples}. The default value is \code{TRUE} if any of
##'     the vectors \code{loc}, \code{scale} and \code{shape} has
##'     length \code{> 1} and \code{FALSE} otherwise.
##'
##' 
##' @param x,q Vector of quantiles.
##'
##' @param p Vector of probabilities.
##'
##' @param n Sample size.
##' 
##' @param lower.tail Logical; if \code{TRUE} (default), probabilities
##'     are \eqn{\textrm{Pr}[X \leq x]}{Pr[X <= x]}, otherwise,
##'     \eqn{\textrm{Pr}[X>x]}{Pr[X > x]}.
##'
##' @return A numeric vector with length \code{n} as described in the
##'     \bold{Details} section. When \code{deriv} is \code{TRUE}, the
##'     returned value has an attribute named \code{"gradient"} which
##'     is a matrix with \eqn{n} lines and \eqn{3} columns containing
##'     the derivatives. A row contains the partial derivatives of the
##'     corresponding element w.r.t. the three parameters \code{loc}
##'     \code{scale} and \code{shape} in that order.
##'
##' @details Each of the probability function normally requires two
##'     formulas: one for the non-zero shape case \eqn{\xi \neq 0}{\xi
##'     != 0} and one for the zero-shape case \eqn{\xi = 0}. However
##'     the non-zero shape formulas lead to numerical instabilities
##'     near \eqn{\xi = 0}, especially for the derivatives
##'     w.r.t. \eqn{\xi}. This can create problem in optimization
##'     tasks. To avoid this, a Taylor expansion w.r.t. \eqn{\xi} is
##'     used for \eqn{|\xi| < \epsilon} for a small positive
##'     \eqn{\epsilon}.  The expansion has order \eqn{2} for the
##'     functions (log-density, distribution and quantile), order
##'     \eqn{1} for their first-order derivatives and order \eqn{0}
##'     for the second-order derivatives.
##'
##'     For the \code{d}, \code{p} and \code{q} functions, the GEV
##'     parameter arguments \code{loc}, \code{scale} and \code{shape}
##'     are recycled in the same fashion as the classical R
##'     distribution functions in the \pkg{stats} package, see e.g.,
##'     \code{\link[stats]{Normal}}, \code{\link[stats]{GammaDist}}, ...
##'     Let \code{n} be the maximum length of the four arguments:
##'     \code{x} \code{q} or \code{p} and the GEV parameter arguments,
##'     then the four provided vectors are recycled in order to have
##'     length \code{n}. The returned vector has length \code{n} and
##'     the attributes \code{"gradient"} and \code{"hessian"}, when
##'     computed, are arrays wich dimension: \code{c(1, 3)} and
##'     \code{c(1, 3, 3)}.
##' 
##' @examples
##' ti <- 1:10; names(ti) <- 2000 + ti
##' mu <- 1.0 + 0.1 * ti
##' ## simulate 40 paths
##' y <- rGEV(n = 40, loc = mu, scale = 1, shape = 0.05)
##' matplot(ti, y, type = "l", col = "gray")
##' lines(ti, apply(y, 1, mean))
dGEV <- function(x, loc = 0.0, scale = 1.0, shape = 0.0, log = FALSE,
                 deriv = FALSE, hessian = FALSE) {
    
    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'deriv' is equal to TRUE")
    }

    if (!all(is.finite(loc)) || !all(is.finite(scale)) || !all(is.finite(shape))) {
        stop("GEV parameters must be finite (non NA)")
    }
    
    res <- .Call(Call_dGEV,
                 as.double(x),
                 as.double(loc),
                 as.double(scale),
                 as.double(shape),
                 as.integer(log),
                 as.integer(deriv),
                 as.integer(hessian))
    
    names(res) <- names(x)
    
    if (deriv) {
        
        n <- length(res)
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 3L),
                  dimnames = list(names(x), c("loc", "scale", "shape")))
        
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 3L, 3L),
                      dimnames = list(names(x),
                                      c("loc", "scale", "shape"),
                                      c("loc", "scale", "shape")))
        }
        
    }

    res
    
}

##' @rdname GEV
##' @export
pGEV <- function(q, loc = 0, scale = 1, shape = 0, lower.tail = TRUE,
                 deriv = FALSE) {

    if (!all(is.finite(loc)) || !all(is.finite(scale)) || !all(is.finite(shape))) {
        stop("GEV parameters must be finite (non NA)")
    }

    res <- .Call(Call_pGEV,
                 as.double(q),
                 as.double(loc),
                 as.double(scale),
                 as.double(shape),
                 as.integer(lower.tail),
                 as.integer(deriv))
    
    names(res) <- names(q)
    
    if (deriv) {
        n <- length(res)
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 3),
                  dimnames = list(names(q), c("loc", "scale", "shape")))
    }
    
    res
    
}

##' @rdname GEV
##' @export
qGEV <- function(p, loc = 0.0, scale = 1.0, shape = 0.0, lower.tail = TRUE,
                 deriv = FALSE, hessian = FALSE) {

    if (!all(is.finite(loc)) || !all(is.finite(scale)) || !all(is.finite(shape))) {
        stop("GEV parameters must be finite (non NA)")
    }
    if (min(p, na.rm = TRUE) < 0.0 || max(p, na.rm = TRUE) > 1.0) 
        stop("`p' must contain probabilities in [0, 1]")
 
    if (hessian && !deriv) {
        stop("'hessian' can be TRUE only when 'deriv' is equal to TRUE")
    }
    
    res <- .Call(Call_qGEV,
                 as.double(p),
                 as.double(loc),
                 as.double(scale),
                 as.double(shape),
                 as.integer(lower.tail),
                 as.integer(deriv),
                 as.integer(hessian))
    
    names(res) <- names(p)
    
    if (deriv) {
        n <- length(res)
        attr(res, "gradient") <-
            array(attr(res, "gradient"),
                  dim = c(n, 3),
                  dimnames = list(names(p), c("loc", "scale", "shape")))
        
        if (hessian) {
            attr(res, "hessian") <-
                array(attr(res, "hessian"),
                      dim = c(n, 3L, 3L),
                      dimnames = list(names(p),
                                      c("loc", "scale", "shape"),
                                      c("loc", "scale", "shape")))
        }
    }
    
    res
    
}

##' @rdname GEV
##' @export
rGEV <- function(n, loc = 0.0, scale = 1.0, shape = 0.0, array) {

    if (!all(is.finite(loc)) || !all(is.finite(scale)) || !all(is.finite(shape))) {
        stop("GEV parameters must be finite (non NA)")
    }
    
    if (missing(array)) {
        array <- max(c(length(loc), length(scale), length(shape))) > 1L
    }
    
    if (array) {

        ## If 'array' is TRUE we return a matrix having the simulated values
        ## as its columns
        L <- .reshapeGEV(x = 1.0, loc = loc, scale = scale, shape = shape)
        
        res <- array(NA, dim = c(L$n, n),
                   dimnames = list(L$nms, paste("sim", 1:n, sep = "")) )
        
        loc <- array(L$loc, dim = c(L$n, n))
        scale <- array(L$scale, dim = c(L$n, n))
        shape <- array(L$shape, dim = c(L$n, n))

        nn <- L$n * n
        res <- array(qGEV(runif(nn), loc = loc, scale = scale, shape = shape),
                     dim = c(L$n, n),
                     dimnames = list(NULL, paste("sim", 1:n, sep = "")))
        
    } else {
        res <- qGEV(runif(n), loc = loc, scale = scale, shape = shape)
    }

    res 
    
}

