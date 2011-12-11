\name{quantileCut}
\alias{quantileCut}

\title{
Cut by quantiles
}

\description{Cuts a variable into equal sized categories
}

\usage{
quantileCut(x,n,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A vector containing the observations.}
  \item{n}{Number of categories}
  \item{...}{Additional arguments to cut}
}

\details{ It is sometimes convenient (though not always wise) to split a continuous numeric variable \code{x} into a set of \code{n} discrete categories that contain an approximately equal number of cases. The \code{quantileCut} function does exactly this. The actual categorisation is done by the \code{\link{cut}} function. However, instead of selecting ranges of equal sizes (the default behaviour in \code{cut}), the \code{quantileCut} function uses the \code{\link{quantile}} function to select  unequal sized ranges so as to ensure that each of the categories contains the same number of observations. The intended purpose of the function is to assist in exploratory data analysis; it is not generally a good idea to use the output of \code{quantileCut} function as a factor in an analysis of variance, for instance, since the factor levels are not interpretable and will almost certainly violate homogeneity of variance. 
}

\value{A factor containing \code{n} levels. The factor levels are determined in the same way as for the \code{cut} function, and can be specified manually using the \code{labels} argument, which is passed to the \code{cut} function.}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Daniel Navarro
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cut}}, \code{\link{quantile}}
}


\examples{
### An example illustrating why care is needed ###

dataset <- c( 0,1,2, 3,4,5, 7,10,15 )       # note the uneven spread of data
x <- quantileCut( dataset, 3 )              # cut into 3 equally frequent bins
table(x)                                    # tabulate
#
# (-0.015,2.67]   (2.67,5.67]     (5.67,15] 
#             3             3             3
# 

# Notice the uneven bin sizes: category 1 covers a range from 0 to 2.67 and 
# category 2 covers a similarly sized range from 2.67 to 5.67, but the third 
# category covers a much larger range, from 5.67 to 15. These categories might
# be useful in some contexts (e.g., the data are ordinal scale), but it is
# important to check that this is so.

# For comparison purposes, here is the behaviour of the more standard cut 
# function when applied to the same data:

y <- cut( dataset, 3 )
table(y)
#
# (-0.015,5]     (5,10]    (10,15] 
#          5          3          1
#

# This time  the categories cover an equal range but have highly unequal
# frequencies.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
