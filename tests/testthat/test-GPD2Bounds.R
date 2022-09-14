context("GPD2Bounds")

## ***************************************************************************
## AUTHOR: Yves Deville <deville.yves@alpestat.com>
## 
## GOAL: Test that the end-points of the GPD2 distribution (finite or
## infinite) are OK
##
## ***************************************************************************

library(numDeriv)
library(testthat)

set.seed(1234)

xi <- rnorm(1, mean = 0, sd = 0.1)
sigma <- rlnorm(1, meanlog = 0, sdlog = 1)

deriv <- TRUE
hessian <- TRUE

lowerEP <- 0.0
if (xi > 0) {
    upperEP <- Inf
    xout <- seq(from = lowerEP - 3 * sigma, to = lowerEP, length = 100)
} else if (xi < 0) {
    upperEP <- - sigma / xi
    xout <- seq(from = upperEP, to =  upperEP + 3 * sigma, length = 100)    
}

quant0 <-  qGPD2(0.0, scale = sigma, shape = xi,
                deriv = deriv, hessian = hessian)

test_that(desc = sprintf("lower end-point"),
          expect_true(quant0 == lowerEP))

quant1 <- qGPD2(1.0, scale = sigma, shape = xi,
               deriv = deriv, hessian = hessian)

test_that(desc = sprintf("upper end-point"),
          expect_true(quant1 == upperEP))

test_that(desc = sprintf("distribution function at -Inf"),
          expect_true(pGPD2(-Inf, scale = sigma, shape = xi,
                           deriv = deriv) == 0.0))

test_that(desc = sprintf("distribution function a -Inf"),
          expect_true(pGPD2(Inf, scale = sigma, shape = xi,
                           deriv = deriv) == 1.0))

test_that(desc = sprintf("distribution function at lower end-point"),
          expect_true(pGPD2(lowerEP, scale = sigma, shape = xi,
                           deriv = deriv) == 0.0))

test_that(desc = sprintf("distribution function at upper end-point"),
          expect_true(pGPD2(upperEP, scale = sigma, shape = xi,
                            deriv = deriv) == 1.0))

densOut <- dGPD2(xout, scale = sigma, shape = xi,
                                     deriv = deriv, hessian = hessian)

test_that(desc = sprintf("density outside of the support"),
          expect_true(all(densOut < 1e-16)))

test_that(desc = sprintf("gradient of the density outside of the support"),                      
          expect_true(all(abs(attr(densOut, "gradient")) < 1e-16)))

test_that(desc = sprintf("Hessian of the density outside of the support"),  
          expect_true(all(abs(attr(densOut, "hessian")) < 1e-16)))
          
