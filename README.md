![R-CMD-check](https://github.com/yvesdeville/nieve/workflows/actions/R-cmd-check/badge.svg)

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

-   Provides probability functions which are *vectorised w.r.t. the
    parameters*, as required in non-stationary EV models or in Bayesian
    inference.

-   Provides the *transformations for the two usual parameterisations*
    of Peaks Over Threshold (POT) models: *Poisson-GP* and *Point
    Process* (PP). The two transformations come with their derivatives,
    as needed to compute the covariance matrices.

Although several R packages devoted EVA compute the exact derivatives
w.r.t. the parameters (**extRemes**, **mev**, …), to our best knowledge,
none of these packages make the derivatives available via exported
functions. See the vignette [**nieve**: Yet Another Extreme Value
Package?](https://github.com/yvesdeville/nieve/blob/main/vignettes/nieve.pdf)
for more details.

Note that **nieve** has no special relation with the famous package
**snow**:)

Install release version from GitHub
===================================

Using the *devtools* package
----------------------------

Note that if you are using Windows, you need to have the
[Rtools](https://cran.r-project.org/bin/windows/Rtools) installed.
Provided that the **devtools** package is installed you can then in an R
session use

    library(devtools)
    install_github("yvesdeville/nieve", dependencies = TRUE, auth_token = myToken)

where `myToken` stands for *your* token. This should install the package
and make it ready to use.

You can also select a specific branch or a specific commit by using the
suitable syntax for `install_github`, see the **devtools** package
documentation.

Clone, build and install
------------------------

### Cloning the repository

If you do not have yet a local `nieve` repository, use `git clone` to
clone the `nieve` repository

    git clone https://github.com/yvesdeville/nieve

This will create a `nieve` sub-directory of the current directory,
i.e. the directory from which the git command was issued. Of course this
can work only if you have the authorisation to clone.

### Installation on Unix and MacOs systems

With these systems you can install a package from its source. Move to
the parent directory of your cloned repository and use the following
command from a terminal to create a tarball source file

    R CMD build nieve

This will produce a source tarball `nieve_x.y.z` where `x`, `y` and `z`
stand for the major, minor and patch version numbers. Then you can
install from a command line

    R CMD INSTALL nieve_x.y.z

Note that you must also have all the packages required by **nieve**
installed.

If you are using the **RStudio** IDE, you can alternatively use menus.

### Install and precompile for Windows

In order to install the package from its source, you must have a
suitable Windows platform with
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) installed. Then
you can proceed as Unix or MacOS users, with a `build` step from command
line.

If you can not (or do not want to) install the **Rtools** you may get a
trusted binary from a friend or colleague next to you.

### Creating a binary precompiled file for Windows

If you have the **Rtools** installed, you can create a binary. Using a
terminal, move if necessary by using `cd` to the directory containing
the source tarball and R command, and then type

    R CMD INSTALL --build nieve_x.y.z

This will create a `.zip` file that can be used on a Windows platform
which may not be equipped with *Rtools*. For instance, with **RStudio**
you can use the menu `Tools/Install Packages` and select
`Install from:`.

### Precompiled binaries

You can make the resulting binary file available to Windows users who do
not have the **Rtools** installed and are allowed to use the package.

Make sure that the **major** version number `x` and the minor version
number `y` are the same as those of the R environment where the package
is to be used.
