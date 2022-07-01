TODO list for the **yaev** package
================
Yves Deville <deville.yves@alpestat.com>

## Package management

-   Use the
    [rchk](https://developer.r-project.org/Blog/public/2019/04/18/common-protect-errors/)
    tool to check for possible memory errors/problems before submitting
    to the CRAN.

-   Make **NSGEV** and **potomax** rely on **yaev** for the distribution
    functions and discard the implementation therein.

-   Make sure that the the `inst/doc` directory is not provided to the
    CRAN since it may fail to conform to CRAN rules. Nearly 1MO of pdf
    files:( Suitable links should be found in the package manual. Maybe
    keep `yaev.Rnw` as a package vignette?

## Computing

-   Implement the Taylor approximation for
    ![\\xi \\approx 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi%20%5Capprox%200 "\xi \approx 0")
    in the `PP2poisGP`and `poisGP2PP` (C) functions.

-   Add the Extrended GPD distributions of Papastathopoulos & Tawn.

## Testing and polishing

-   Remove any kind of dependence to **Renext**.

-   Should we keep the R implementation for the GEV distribution?

-   Add more tests for the case
    ![\\xi \\approx 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cxi%20%5Capprox%200 "\xi \approx 0")
    in the transformations. `PP2poisGP`and `poisGP2PP`.

-   Make the quantile functions work correctly for `p = 0.0` and
    `p =   1.0`. The corresponding end-point of the distribution (be it
    finite or infinite) must be returned.

-   Make sure that the distribution functions work correctly with
    `q = -Inf` and `q = Inf`. The value `0` or `1` must be returned.

-   Make sure that the densities are zero (not `NA`) outside of the
    support.

-   Make sure that the probability functions return `NA` when an
    incorrect parameter is used (negative `scale`).
