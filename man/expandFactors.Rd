\name{expandFactors}
\alias{expandFactors}

\title{Expand factors to a set of contrasts}

\description{
Substitutes all factors in a data frame with the set of contrasts with which that factor is associated
}

\usage{
	expandFactors( data, ... )
}

\arguments{
  \item{data}{A data frame.}
  \item{...}{Additional arguments to be passed to model.matrix}
}

\details{
The \code{expandFactors} function replaces all of the factors in a data frame with the set of contrasts output by the \code{contrasts} function or \code{model.matrix}. It may be useful for teaching purposes when explaining relationship between ANOVA and regression. 
}

\value{
A data frame.
}



%\references{
%% ~put references to the literature/web site here ~
%}

\author{Daniel Navarro}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{model.matrix}}, \code{\link{contrasts}}
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

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
