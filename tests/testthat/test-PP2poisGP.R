context("PP2poisGP")

## ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## GOAL: Test the transformation PP -> Poisson-GP and its derivative
## ***************************************************************************
library(numDeriv)
library(testthat)

set.seed(1234)

n <- 10
locStar <- rnorm(n, sd = 20)
scaleStar <- rgamma(n, shape = 2)
shapeStar <- rnorm(n, sd = 0.1)

p <- runif(n)
threshold <- rep(NA, n)
for (i in 1:n) {
    threshold[i] <- qGEV(p = p[i], loc = locStar[i], scale = scaleStar[i],
                         shape = shapeStar[i])
}

## ===========================================================================
## Check the n cases with non-zero shape
## ===========================================================================

for (i in n) {
    
    res <- PP2poisGP(locStar = locStar[i],
                     scaleStar = scaleStar[i],
                     shapeStar = shapeStar[i],
                     threshold = threshold[i],
                     deriv = TRUE)
    
    thetaStar <- c("loc" = locStar[i], "scale" = scaleStar[i],
                   "shape" = shapeStar[i])
    
    ## Compute with Renext
    res0 <- Renext::gev2Ren(thetaStar, threshold = threshold[i],
                            distname.y = "GPD")
    
    g <- drop(attr(res, "gradient"))
    g0 <- attr(res0, "jacobian")
    
    attr(res, "gradient") <- NULL
    attr(res0, "jacobian") <- attr(res0, "threshold") <- NULL
    
    e <- drop(res) - res0
    test_that(desc = sprintf("Case shape %s, consistency of value with Renext",
                  "non-zero"),
              expect_lt(max(abs(e)), 1e-10))

    e <- g - g0[-2, ]

    test_that(desc = sprintf("Case shape %s, consistency of jacobian with Renext",
                  "non-zero"),
              expect_lt(max(abs(e)), 1e-10))

}

## ===========================================================================
## Check case 'n' zero shape. Use the Jacobian
## ===========================================================================

shapeStar <- runif(n, min = -1e-5, max = 1e-5)

PP2poisGPFun <- function(thetaStar, threshold) {
    PP2poisGP(locStar = thetaStar[1],
              scaleStar = thetaStar[2],
              shapeStar = thetaStar[3],
              threshold = threshold,
              deriv = FALSE)
    
}

for (i in 1:n) {

    thetaStar <- c("loc" = locStar[i], "scale" = scaleStar[i],
                   "shape" = shapeStar[i])
    
    jacNum <- jacobian(func = PP2poisGPFun, x = thetaStar, threshold = threshold[i])
    
    res <- PP2poisGP(locStar = locStar[i],
                 scaleStar = scaleStar[i],
                 shapeStar = shapeStar[i],
                 threshold = threshold[i],
                 deriv = TRUE)
    g <- drop(attr(res, "gradient"))
    
    e <- g - jacNum
    
    cond <- (max(abs(g - jacNum)) < 4e-2) ||
        (max(abs(g - jacNum) / (abs(g) + 1e-9)) < 4e-2)
    
    if (!cond) {
        print(g)
        print(jacNum)
    }
    
    test_that(desc = sprintf("Case shape %s, Jacobian and numeric Jacobian",
                             "zero"),
              expect_true(cond))

}
