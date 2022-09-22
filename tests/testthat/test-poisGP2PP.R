context("poisGP2PP")

## ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## GOAL: Test the transformation poisGP -> PP and its derivative
## ***************************************************************************

library(numDeriv)
library(testthat)

set.seed(1234)

## =============================================================================
## We begin with 'n' checks for the case where the shape is
## non-zero. We then check the consistency with Renext
## =============================================================================

n <- 6
lambda <- rexp(n)
loc <- rnorm(n, sd = 20)
scale <- rgamma(n, shape = 2)
shape <- rnorm(n, sd = 0.1)
w <- rexp(n)

for (i in seq_along(shape)) {
    
    res <- poisGP2PP(lambda = lambda[i],
                     loc = loc[i], scale = scale[i], shape = shape[i],
                     deriv = TRUE)
    
    theta <- c("lambda" = lambda[i], "scale" = scale[i], "shape" = shape[i])
    
    ## Compute with Renext
    res0 <- Renext::Ren2gev(theta, threshold = loc[i])
    
    g <- drop(attr(res, "gradient"))
    g0 <- attr(res0, "jacobian")
    
    attr(res, "gradient") <- NULL
    attr(res0, "jacobian") <- attr(res0, "threshold") <- NULL
    
    e <- drop(res) - res0
    test_that(desc = sprintf("Case shape %s, consistency of value with Renext",
                  "non-zero"),
              expect_lt(max(abs(e)), 1e-10))
    
    e <- g[ , -2]  - g0

    test_that(desc = sprintf("Case shape %s, consistency of jacobian with Renext",
                  "non-zero"),
              expect_lt(max(abs(e)), 1e-10))

}

## ===========================================================================
## Check case 'n' zero shape. Use the Jacobian
##
## Note that it was found that the numeric derivatives suffer from a
## lack of precision (difficult to explain). So we decided to compare the
## approximation and the true formula as implemented in 'Renext'.
## ===========================================================================

shape <- runif(n, min = -1e-5, max = 1e-5)

poisGP2PPFun <- function(theta, loc) {
    
    poisGP2PP(lambda = theta[1],
              scale = theta[2],
              shape = theta[3],
              loc = loc,
              deriv = FALSE)
    
}

for (i in 1:n) {

    theta <- c("lambda" = lambda[i], "scale" = scale[i],
               "shape" = shape[i])
    
    jacNum <- jacobian(func = poisGP2PPFun, x = theta, loc = loc[i])
    
    res <- poisGP2PP(lambda = lambda[i],
                     loc = loc[i],
                     scale = scale[i],
                     shape = shape[i],
                     deriv = TRUE)
    g <- drop(attr(res, "gradient"))
    g <- g[ , -2]
    
    res0 <- Renext::Ren2gev(theta, threshold = loc[i],
                            distname.y = "GPD")
    
    jacTest <- attr(res0, "jacobian")
    
    cond <- (max(abs(g - jacTest)) < 1e-3) ||
        (max(abs(g - jacTest) / (abs(g) + 1e-9)) < 1e-3)
    
    if (!cond) {
        cat("================= test failure details ================\n")
        cat("o Parameter value 'theta'\n")
        print(theta)
        attributes(res) <- NULL
        attributes(res0) <- NULL
        fval <- rbind("nieve" = res, "Renext" = res0)
        colnames(fval) <- c("locStar", "scaleStar", "shapeStar")
        cat("o Function value\n")
        print(fval)
        cat("\no Jacobian matrix\n")
        colnames(jacTest) <- paste0(c("loc", "scale", "shape"), "Test")
        print(cbind(g, jacTest))
        cat("\n")
    }
    
    test_that(desc = sprintf("Case shape %s, Jacobian and numeric Jacobian",
                             "zero"),
              expect_true(cond))

}
