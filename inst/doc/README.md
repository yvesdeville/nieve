Compile a LaTeX document using Maxima via the **maxiplot** package
------------------------------------------------------------------

In order to compile a LaTeX doucment with maxima

    pdflatex mydocument.tex

    maxima -b mydocument.mac

Compiling whole document
------------------------

The documentations comes as 4 documents. The document `N_` do not have
any header and can not compiled

-   `GEV.tex` GEV distribution functions and their derivatives.
-   `GPD2.tex` the two-parameter GP2 distribution functions and their
    derivatives.
-   `PP2PoisGP.tex` The PP to Poisson-GP transform and its derivatives.
-   `PoisGP2PP.tex` The Poisson-GP to PP transform and its derivatives.
