<!-- badges: start -->

[![R-cmd-check](https://github.com/yvesdeville/nieve/actions/workflows/check.yml/badge.svg)](https://github.com/yvesdeville/nieve/actions/workflows/check.yml)
[![Codecov test
coverage](https://codecov.io/gh/yvesdeville/nieve/branch/main/graph/badge.svg)](https://app.codecov.io/gh/yvesdeville/nieve?branch=main)
[![CRAN
status](https://www.r-pkg.org/badges/version/nieve)](https://cran.r-project.org/package=nieve)
[![Downloads
(monthly)](https://cranlogs.r-pkg.org/badges/nieve?color=brightgreen)](https://cran.r-project.org/package=nieve)
[![Downloads
(total)](https://cranlogs.r-pkg.org/badges/grand-total/nieve?color=brightgreen)](https://cran.r-project.org/package=nieve)
<!-- badges: end -->

<img src="inst/images/nieve.png" height = "150" align="center"/>

Goals and scope
===============

The **nieve** package was partly funded by the French [*Institut de
Radioprotection et Sûreté Nucléaire* (IRSN)](https://www.irsn.fr/) and
some of the code formerly was part of R packages owned by IRSN/Behrig.

The **nieve** package is intended to be a “low-level” package, providing
fast and well-tested “basic” functions for Extreme Value Analysis (EVA).
It is not intended to provide sophisticated EVA models which should be
found or be implemented in other packages.

The package

-   Provides the *derivatives w.r.t. the parameters* for the probability
    functions related to the Generalized Pareto (GP) and the Generalized
    Extreme Value (GEV) distributions. This involves: the log-density,
    the distribution function or survival and the quantile function. The
    2-nd order derivative (Hessian) is available in most cases.

-   Provides probability functions which are *vectorized w.r.t. the
    parameters*, as required in non-stationary EV models or in Bayesian
    inference.

-   Provides the *transformations for the two usual parameterizations*
    of Peaks Over Threshold (POT) models: *Poisson-GP* and *Point
    Process* (PP). The two transformations come with their derivatives,
    as needed to compute the covariance matrices.

Although several R packages devoted EVA compute the exact derivatives
w.r.t. the parameters (**extRemes**, **mev**, …), to our best knowledge,
none of these packages make the derivatives available via exported
functions. See the vignette [**nieve**: Yet Another Extreme Value
Package?](https://github.com/yvesdeville/nieve/blob/main/vignettes/nieve.pdf)
for more details.

The vignette [Special Values in Extreme-Value
Distributions](https://github.com/yvesdeville/nieve/blob/main/inst/noncranVignettes/nonFiniteValues.pdf)
compares the treatment of non-finite values in the probability functions
across Extreme-Value packages.

Note that **nieve** has no special relation with the famous package
**snow**:)

News
====

See the
[NEWS.md](https://github.com/yvesdeville/nieve/blob/main/NEWS.md) file
or (after install) use `news(package = "nieve")`.

Installation
============

### From CRAN

**nieve** is on CRAN since 2022-12-02 hence can be installed using the
`install.packages` function or the **RStudio** IDE.

### From GitHub

Note that this package needs a compilation since it embeds some code
written in C. If you are using Windows, you need to have the
[Rtools](https://cran.r-project.org/bin/windows/Rtools) installed.

Provided that the **devtools** package is installed you can then in an R
session use

    devtools::install_github("yvesdeville/nieve", dependencies = TRUE)

You can also select a specific branch or a specific commit by using the
suitable syntax for `install_github`, see the **devtools** package
documentation.

Alternatively you can Clone the repository using `git clone`, and then
build and install.
