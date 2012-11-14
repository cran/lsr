\name{xfun}

\alias{xfun}

\title{Apply a function by multiple grouping variables}

\description{
Apply a function by multiple grouping variables
}

\usage{
	xfun( formula, data = NULL, fun, subset = NULL, split = TRUE, ... )
}

\arguments{
  \item{formula}{A two-sided formula specifying the outcome variable and the grouping variables}
  \item{data}{An optional data frame containing the variables}
  \item{fun}{The function to be applied to the outcome}
  \item{subset}{A vector used to specify a subset of the cases}
  \item{split}{Should \code{xfun} attempt to "split" an array of lists into a list of arrays? (See Details)}
  \item{...}{Additional arguments to be passed to \code{tapply} or to the functions specified by \code{fun}}
}

\details{
The \code{xfun} function applies the function \code{fun} to an outcome variable, broken down by one or more grouping variables. The function is a wrapper to \code{tapply}, but the usage is similar to the \code{xtabs} function used to construct cross-tabulations. The outcome variable and grouping variables are specified by a two-sided \code{formula} of the form \code{outcome ~ group1 + group2 + ... + groupK}, where the variables are stored in the data frame specified by \code{data}. The optional \code{subset} argument can be used to subset the data frame if need be (see Example 2). The function to be applied to the outcome variable is specified by the \code{fun} argument. 

The simplest way to use \code{xfun} is to compute the mean, standard deviation or some other function that returns a single value separately for each group (see Example 1). If there K grouping variables, the output in such cases will be a K-dimensional array, such that each element in the array contains the output of \code{fun} when applied to the data for the corresponding group. 

If the function \code{fun} returns a more complex object, then the behaviour of \code{xfun} is more complicated. If \code{split=TRUE} (the default) and if the output of \code{fun} is always a vector, then \code{xfun} will return a list, such that the first element of the list is an array containing the first element of each vector output by \code{fun}, the second element of the list is an array containing the second element of each vector, etc. For example, if \code{fun=ciMean}, then the output is a list of length 2: the first element of the list is an array containing all the lower bounds of the confidence intervals, and the second element of the list is an array containing all the upper bounds of the confidence intervals (see Examples 3 and 5). If \code{split=FALSE} or if the output of \code{fun} is not a vector for at least one group, then \code{xfun} reverts to the behaviour of \code{tapply}, and produces an array of mode list (see Example 4).
}

\value{
Depending on what the function \code{fun} returns and on the value of the \code{split} argument, the output is either an array or a list (see Details).
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
 \code{\link{tapply}},\code{\link{xtabs}}
}

\examples{

#### Example 1: basic usage

# data
df <- data.frame( outcome = 1:6, group1 = c(1,1,1,2,2,2), group2 = c(1,2,1,2,1,2) )

# use xtabs() to obtain a frequency table for group1 x group2:
xtabs( ~ group1 + group2, df )
#       group2
# group1 1 2
#      1 2 1
#      2 1 2

# use xfun() to obtain the same frequency table:
 xfun( outcome ~ group1 + group2, df, length )
#       group2
# group1 1 2
#      1 2 1
#      2 1 2

# use xfun() to find the group means:
xfun( outcome ~ group1 + group2, df, mean )
xfun( formula = outcome ~ group1 + group2, data = df, fun = mean )
#       group2
# group1 1 2
#      1 2 2
#      2 5 5

# use xfun() to find the smallest value in each group:
xfun( outcome ~ group1 + group2, df, min )
#       group2
# group1 1 2
#      1 1 2
#      2 5 4


#### Example 2: subsetting the data frame

xfun( formula = outcome ~ group1 + group2, data = df, fun = min, subset = -(5:6) )
#       group2
# group1  1 2
#      1  1 2
#      2 NA 4


#### Example 3: by default, xfun produces a list of arrays 
#### when fun produces a vector of outputs 

df2 <- data.frame( outcome = 1:12, group1 = gl(3,4), group2 = gl(2,1,12) )
xfun( outcome ~ group1 + group2, df2, range )
xfun( outcome ~ group1 + group2, df2, ciMean )


#### Example 4: if split=FALSE, or if the output of fun is 
#### not a vector, then xfun reverts to an array of mode list, 
#### as per the tapply function:

xfun( outcome ~ group1 + group2, df2, ciMean, split=FALSE )


#### Example 5: if fun returns outputs of different lengths 
#### for different groups, the remaining entries are padded 
#### with NAs:

xfun( formula = outcome ~ group1 + group2, data = df, fun = c )

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
