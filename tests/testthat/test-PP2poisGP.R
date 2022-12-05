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
shapeStar <- rnorm(n, sd = 5e-3)

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
    cond <- max(abs(e)) < 1e-8 
    if (!cond) {
        cat("\nFunction value: nieve and Renext\n")
        print(shapeStar[i])
        print(res)
        print(res0)
    }
    test_that(desc = sprintf("Case shape %s, consistency of value with Renext",
                  "non-zero"),
              expect_true(cond))

    e <- g - g0[-2, ]
    cond <- max(abs(e)) < 1e-8
    if (!cond) {
        cat("\nGradient: nieve and Renext\n")
        print(shapeStar[i])
        print(g)
        print(g0[-2, ])
    }
    test_that(desc = sprintf("Case shape %s, consistency of jacobian with Renext",
                  "non-zero"),
              expect_true(cond))

}

## ===========================================================================
## Check case 'n' zero shape. Use the Jacobian
##
## Note that it was found that the numeric derivatives suffer from a
## lack of precision (difficult to explain). So we decided to compare the
## approximation and the true formula as implemented in 'Renext'.
## ===========================================================================

shapeStar <- runif(n, min = -1e-4, max = 1e-4)

PP2poisGPFun <- function(thetaStar, threshold) {
    res <- PP2poisGP(locStar = thetaStar[1],
                     scaleStar = thetaStar[2],
                     shapeStar = thetaStar[3],
                     threshold = threshold,
                     deriv = FALSE)
}


for (i in 1:n) {

    thetaStar <- c("loc" = locStar[i], "scale" = scaleStar[i],
                   "shape" = shapeStar[i])

    jacNum <- jacobian(func = PP2poisGPFun, x = thetaStar,
                       threshold = threshold[i],
                       method = "simple", method.args = list(eps = 6e-3))

    res0 <- Renext::gev2Ren(thetaStar, threshold = threshold[i],
                            distname.y = "GPD")
    
    jacTest <- attr(res0, "jacobian")[-2, ]
    
    ## method = "simple", method.args = list(eps = 1e-2))
 
    res <- PP2poisGP(locStar = locStar[i],
                     scaleStar = scaleStar[i],
                     shapeStar = shapeStar[i],
                     threshold = threshold[i],
                     deriv = TRUE)
    
    g <- drop(attr(res, "gradient"))
    
    cond <- (max(abs(g - jacTest)) < 1e-3 / PREC) ||
        (max(abs(g - jacTest) / (abs(g) + 1e-9)) < 1e-3 / PREC)
    
    if (!cond) {
        cat("================= test failure details ================\n")
        cat("o Parameter value 'thetaStar'\n")
        print(thetaStar)
        cat("\no Jacobian matrix\n")
        colnames(jacNum) <- paste0(c("loc", "scale", "shape"), "Test")
        print(cbind(g, jacTest))
        cat("\n")
    }
    
    test_that(desc = sprintf("Case shape %s, Jacobian and numeric Jacobian",
                             "zero"),
              expect_true(cond))

}


