TODO list for the **yaev** package
----------------------------------

-   Implement the Taylor approximation for *ξ* ≈ 0 in the `PP2poisGP`and
    `poisGP2PP` (C) functions.

-   Use the
    [rchk](https://developer.r-project.org/Blog/public/2019/04/18/common-protect-errors/)
    tool to check for possible memory erros before submitting to the
    CRAN.

-   Make **NSGEV** and **potomax** rely on **yaev** for the distribution
    functions.

-   Add more tests for the case *ξ* ≈ 0 in the transformations.
    `PP2poisGP`and `poisGP2PP`.

-   Add the Extrended GPD distributions of Papastathopoulos & Tawn.
