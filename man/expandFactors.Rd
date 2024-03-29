% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expandFactors.R
\name{expandFactors}
\alias{expandFactors}
\title{Expand factors to a set of contrasts}
\usage{
expandFactors(data, ...)
}
\arguments{
\item{data}{A data frame.}

\item{...}{Additional arguments to be passed to model.matrix}
}
\value{
A data frame.
}
\description{
Substitutes all factors in a data frame with the set of
contrasts with which that factor is associated
}
\details{
The \code{expandFactors} function replaces all of the factors
in a data frame with the set of contrasts output by the \code{contrasts}
function or \code{model.matrix}. It may be useful for teaching purposes
when explaining relationship between ANOVA and regression.
}
\examples{
grading <- data.frame( teacher = factor( c("Amy","Amy","Ben","Ben","Cat") ),
                       gender = factor( c("male","female","female","male","male") ),
                       grade = c(75,80,45,50,65) )

# expand factors using the default contrasts (usually treatment contrasts)
expandFactors( grading )

# specify the contrasts using the contrasts.arg argument to model.matrix
my.contrasts <- list( teacher = "contr.helmert", gender = "contr.treatment" )
expandFactors( grading, contrasts.arg = my.contrasts )

}
