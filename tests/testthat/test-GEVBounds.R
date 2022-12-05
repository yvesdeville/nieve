context("GEVBounds")

## ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## 
## GOAL: Test that the end-points of the GEV distribution (finite or
## infinite) are OK
##
## ***************************************************************************

library(numDeriv)
library(testthat)

set.seed(1234)

xi <- rnorm(1, mean = 0, sd = 0.1)
mu <- rnorm(1, mean = 0, sd = 3)
sigma <- rlnorm(1, meanlog = 0, sdlog = 1)

deriv <- TRUE
hessian <- TRUE

if (xi > 0) {
    lowerEP <- mu - sigma / xi
    upperEP <- Inf
    xout <- seq(from = lowerEP - 3 * sigma, to = lowerEP, length = 100)
} else if (xi < 0) {
    lowerEP <- -Inf
    upperEP <- mu - sigma / xi
    xout <- seq(from = upperEP, to =  upperEP + 3 * sigma, length = 100)
    
}

quant0 <-  qGEV(0.0, loc = mu, scale = sigma, shape = xi,
                deriv = deriv, hessian = hessian)

test_that(desc = sprintf("lower end-point"),
          expect_true(quant0 == lowerEP))

quant1 <- qGEV(1.0, loc = mu, scale = sigma, shape = xi,
               deriv = deriv, hessian = hessian)

test_that(desc = sprintf("upper end-point"),
          expect_true(quant1 == upperEP))

test_that(desc = sprintf("distribution function at -Inf"),
          expect_true(pGEV(-Inf, loc = mu, scale = sigma, shape = xi,
                           deriv = deriv) == 0.0))

test_that(desc = sprintf("distribution function a -Inf"),
          expect_true(pGEV(Inf, loc = mu, scale = sigma, shape = xi,
                           deriv = deriv) == 1.0))

test_that(desc = sprintf("distribution function at lower end-point"),
          expect_true(pGEV(lowerEP, loc = mu, scale = sigma, shape = xi,
                           deriv = deriv) == 0.0))

test_that(desc = sprintf("distribution function at upper end-point"),
          expect_true(pGEV(upperEP, loc = mu, scale = sigma, shape = xi,
                           deriv = deriv) == 1.0))

densOut <- dGEV(xout, loc = mu, scale = sigma, shape = xi,
                deriv = deriv, hessian = hessian)

test_that(desc = sprintf("density outside of the support"),
          expect_true(all(densOut < 1e-16 / PREC)))

test_that(desc = sprintf("gradient of the density outside of the support"),                      
          expect_true(all(abs(attr(densOut, "gradient")) < 1e-16 / PREC)))

test_that(desc = sprintf("Hessian of the density outside of the support"),  
          expect_true(all(abs(attr(densOut, "hessian")) < 1e-16 / PREC)))
          
