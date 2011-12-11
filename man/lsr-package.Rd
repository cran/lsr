\name{lsr-package}
\alias{lsr-package}
\alias{lsr}
\docType{package}

\title{
Companion to "Learning Statistics with R" 
}

\description{
Functions that I have found useful for teaching purposes.
}

\details{
\tabular{ll}{
Package: \tab lsr\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2011-12-11\cr
License: \tab GPL-3 \cr
LazyLoad: \tab yes\cr
}

The lsr package is a companion to the lecture notes from which I teach my introductory statistics class (for psychology students). It takes its name from the name of the lecture notes, "Learning Statistics with R". The purpose of the package is to provide a few functions that are useful for beginners, but not necessarily needed by experienced users. An example of this is the \code{\link{ciMean}} function, which generates confidence intervals around the mean for a normally distributed variable; bypassing the need to introduce the rather more complicated \code{\link{confint}} function early on in the learning process, or to call the \code{\link{t.test}} function before hypothesis testing has been discussed.

Besides the confidence interval function, some functions that I have found particularly handy while teaching are the \code{\link{sortFrame}} function, 
the functions used to "automatically" switch between data frames in wide form versus long form (e.g., \code{\link{wideToLong}}), functions for calculating effect sizes for simple models (e.g., \code{\link{cohensD}}), and the function for calculating the sample mode (\code{\link{modeOf}}).

--- WARNING: This package is still under development. The functions work for most situations that I envisage them being used for, but some extensions are planned. Also note that very little exception handling has been written at this stage. Thus, although the package serves its basic goals (it works for all examples I need for my 2012 statistics class), until this notice disappears some care (that is, more than usual) is advised in the use of this package. ---

}

\author{
Author: Daniel Navarro
\cr
Maintainer: Daniel Navarro <daniel.navarro@adelaide.edu.au>
}

%\references{}

\keyword{ package }

%\seealso{
%}

%\examples{
%## ~~ simple examples of the most important functions ~~
%}
