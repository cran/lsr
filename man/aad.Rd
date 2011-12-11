\name{aad}
\alias{aad}

\title{
Mean (average) absolute deviation from the mean
}

\description{
Calculates the mean absolute deviation from the sample mean. 
}

\usage{
aad(x, na.rm = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector containing the observations.}
  \item{na.rm}{A logical value indicating whether or not missing values should be removed. Defaults to \code{FALSE}}
}

\details{
The \code{aad} function calculates the average (mean) absolute deviation from the mean value of \code{x}, removing \code{NA} values if requested by the user. It exists primarily to simplify the discussion of descriptive statistics during an introductory stats class.
}

\value{
Numeric.
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Daniel Navarro
}

\note{
The more natural acronym is taken by the \code{\link{mad}} function which calculates the median absolute deviation. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{mad}}
}


\examples{
# basic usage
X <- c(1, 3, 6)  # data 
aad(X)           # returns a value of 2

# removing missing data
X <- c(1, 3, NA, 6)   # data 
aad(X)                # returns NA 
aad(X, na.rm = TRUE)  # returns 2 

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
