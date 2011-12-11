\name{cramersV}
\alias{cramersV}

\title{Cramer's V}

\description{Calculate the Cramer's V measure of association
}

\usage{
cramersV(...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{...}{Arguments to be passed to the \code{chisq.test} function.}
}

\details{
Calculates the Cramer's V measure of effect size for a chi-square goodness of fit test. The arguments to the \code{cramersV} function are all passed straight to the \code{chisq.test} function, and should have the same format.
}

\value{A numeric variable with a single element corresponding to the value of V.}

\section{Warning}{This function is under development. Some details may change in later versions.}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Daniel Navarro
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{chisq.test}}, \code{assocstats} (in the \code{vcd} package)
}


\examples{# Consider an experiment with two conditions, each with 100 
# participants. Each participant chooses between one of three
# options. Possible data for this experiment:

condition1 <- c(30, 20, 50) 
condition2 <- c(35, 30, 35)
X <- cbind( condition1, condition2 )
rownames(X) <- c( 'choice1', 'choice2', 'choice3' )
print(X)
#         condition1 condition2
# choice1         30         35
# choice2         20         30
# choice3         50         35
#

# To test the null hypothesis that the distribution of choices
# is identical in the two conditions, we would run a chi-square
# test:

chisq.test(X)
# 	Pearson's Chi-squared test
# 
# data:  X 
# X-squared = 5.0317, df = 2, p-value = 0.0808

# To estimate the effect size we can use Cramer's V:
cramersV( X )  # returns a value of 0.159


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
