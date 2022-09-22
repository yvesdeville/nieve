Package scope and management
----------------------------

-   **\[X\]** Use the
    [rchk](https://developer.r-project.org/Blog/public/2019/04/18/common-protect-errors/)
    tool to check for possible memory errors/problems before submitting
    to the CRAN.

-   Make **NSGEV** and **potomax** rely on **yaev** for the distribution
    functions, and discard the implementation of the probability
    functions therein.

-   Make sure that the the `inst/computing` directory and the four
    vignettes `GEV.Rnw`, `GPD2.Rnw`, `PP2PoisGP.Rnw`, `PoisGP2PP.Rnw`
    are not provided to the CRAN since the package may then fail to
    conform to CRAN rules. Nearly 1MO of pdf files:( Suitable links
    should be found in the package manual. Maybe keep `yaev.Rnw` as a
    package vignette with a few examples?

-   Remove any kind of dependence to **Renext**?

-   Should we keep the R implementation for the GEV distribution?

Documentation
-------------

-   **\[X\]** Inquire about the **extRemes** package and report it in
    the table of the vignette `yaev.Rnw`. The package does not aim at
    exporting probability functions which are rather seen as being for
    internal use.

-   **{X\]** Add a `.bib` file with the citations for the packages and
    cite the packages in the vignette.

-   In `GPD2.tex`, recall how can the derivatives of the distribution
    function *F* are related to those of the cumulated hazard *H*.

Computing
---------

-   **\[X\]** Implement the Taylor approximation for *ξ* ≈ 0 in the
    `PP2poisGP`and `poisGP2PP` (C) functions.

-   Add the *Extended GPD distributions* of Papastathopoulos & Tawn.
    Since the distribution functions or survival are obtained by
    composition of functions, the derivatives should come by chain rule.

-   Should there be a `Gumbel2` distribution?

-   Use a global variable for `eps` in order to share a value across
    functions? Or add a new argument for each function?

Testing and polishing
---------------------

-   **\[X\]** Add tests for the specific case *ξ* ≈ 0 (more precisely,
    *ξ* &lt; *ϵ*) in the transformations `PP2poisGP`and `poisGP2PP`.

-   Add a test checking that the transformation `PP2poisGP`and
    `poisGP2PP` are reciprocal to each other (in both directions).

-   **\[Done except for exp1\]** Make the quantile functions work
    correctly for `p = 0.0` and `p = 1.0`. The corresponding end-point
    of the distribution (be it finite or infinite) must be returned.

-   **\[Done except for exp1\]** Make sure that the distribution
    functions work correctly with `q = -Inf` and `q = Inf`. The value
    `0` or `1` must be returned.

-   Add the formal argument `log.p` to the distribution functions
    `pGPD2` and `pGEV` to better conform to what is done for
    distribution functions in the **stats** package.

-   Add the formal argument `log.p` to the quantile functions `qGPD2`
    and `qGEV` to better conform to what is done for distribution
    functions. This of course implies changing the derivative and the
    Hessian.

-   **\[X\]** Make sure that the density functions are zero (not `NA`)
    outside of the support.

-   **\[X\]** Make sure that the distribution functions are OK outside
    or the support: zero or one.

-   Make sure that the probability functions return `NA` when an
    incorrect parameter is used (negative `scale`).

-   Test the functions related to the `Exp1`distribution.
