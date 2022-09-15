Package management
------------------

-   **\[X\]** Use the
    [rchk](https://developer.r-project.org/Blog/public/2019/04/18/common-protect-errors/)
    tool to check for possible memory errors/problems before submitting
    to the CRAN.

-   Make **NSGEV** and **potomax** rely on **yaev** for the distribution
    functions and discard the implementation therein.

-   Make sure that the the `inst/computing` directory and the four
    vignettes `GEV.Rnw`, `GPD2.Rnw`, `PP2PoisGP.Rnw`, `PoisGP2PP.Rnw`
    are not provided to the CRAN since the package may then fail to
    conform to CRAN rules. Nearly 1MO of pdf files:( Suitable links
    should be found in the package manual. Maybe keep `yaev.Rnw` as a
    package vignette with a few examples?

Documentation
-------------

-   Inquire about the **extRemes** package and report it in the table of
    the vignette `yaev.Rnw`.

-   Add a `.bib` file with the citations for the packages and cite the
    package in the vignette.

Computing
---------

-   Implement the Taylor approximation for *ξ* ≈ 0 in the `PP2poisGP`and
    `poisGP2PP` (C) functions.

-   Add the Extrended GPD distributions of Papastathopoulos & Tawn.

Testing and polishing
---------------------

-   Remove any kind of dependence to **Renext**.

-   Should we keep the R implementation for the GEV distribution?

-   **\[X\]** Add tests for the specific case *ξ* ≈ 0 (more precisely,
    *ξ* &lt; *ϵ*) in the transformations `PP2poisGP`and `poisGP2PP`.

-   Add a test checking that the transformation `PP2poisGP`and
    `poisGP2PP` are reciprocal to each other (in both directions).

-   **\[X except for exp1\]** Make the quantile functions work correctly
    for `p = 0.0` and `p = 1.0`. The corresponding end-point of the
    distribution (be it finite or infinite) must be returned.

-   **\[X execpt for exp1\]** Make sure that the distribution functions
    work correctly with `q = -Inf` and `q = Inf`. The value `0` or `1`
    must be returned.

-   Add the formal argument `log.p` to the distribution functions
    `pGPD2` and `pGEV` to better conform to what is done for
    distribution functions in the **stats** package.

-   Add the formal argument `log.p` to the quantile functions `qGPD2`
    and `qGEV` to better conform to what is done for distribution
    functions. This of course implies changing the derivative and the
    Hessian.

-   **\[Done for `dGPD2`\]** Make sure that the density functions are
    zero (not `NA`) outside of the support.

-   **\[Done for `pGPD2`\]** Make sure that the distribution functions
    are OK outside or the support (zero or one).

-   Make sure that the probability functions return `NA` when an
    incorrect parameter is used (negative `scale`).

-   In `GPD2.tex`, recall how can the derivatives of the distribution
    function *F* are related to those of the cumulated hazard *H*.

-   Test the functions related to the `Exp1`distribution.
