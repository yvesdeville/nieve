Compile a LaTeX document using Maxima via the **maxiplot** package
------------------------------------------------------------------

In order to compile a LaTeX document with maxima you need to have the
**Maxima** system installed (possibly with its sources) and the
**maxiplot** package working with your LaTeX distribution. Then use

    pdflatex mydocument.tex

This create a number of auxiliary files, including a `.mac` file. Then

    maxima -b mydocument.mac

and compile again

    pdflatex mydocument.tex

This should provide the pdf file embeding the results computed by
**Maxima** In case of problems it may be required to clean the auxiliary
files `.aux`, `.mac` and `.mxp`.

Possible Issues
---------------

With my configuration and **Maxima 5.45.1** the line
`load("mactex.lisp")`in the `maxiplot.sty` file generated an error since
the file seems to no longer exist uin **Maxima** source code. I simply
commented this line (near line 70).

Provided documents
------------------

The documentation includes 5 documents along with their LaTeX or Sweave
sources. Note that rhe documents `N_` do not have any header and can not
compiled, and that the `.tex` documents require the **maxiplot**

-   `nieve.pdf` explains why and how the probability distributions are
    evaluated in **nieve**.

-   [`GEV.pdf`](GEV.pdf) provides details about the Generalized Extreme
    Value distribution `GEV`: probability functions and their
    derivatives.

-   [`GPD2.pdf`](GPD.pdf) provide details about the two-parameter
    Genaralized Pareto distribution `GP2`: probability functions and
    their derivatives.

-   [`PP2PoisGP.pdf`](PP2PoisGP) Is about the *PP to Poisson-GP*
    transform and its derivatives.

-   [`PoisGP2PP.pdf`](PoisGP2PP.pdf) Is about the *Poisson-GP to PP*
    transform and its derivatives.
