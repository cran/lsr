\name{cohensD}
\alias{cohensD}

\title{Cohen's d}
\description{Calculates the Cohen's d measure of effect size.}

\usage{ cohensD( x, y = NULL, data = NULL, method = "pooled", mu = 0 ) }
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Either a numeric variable containing the data for group 1, or a formula of the form \code{outcome ~ group}}
  \item{y}{If \code{x} is a numeric variable, the \code{y} argument should be a numeric variable containing the data for group 2. If a one-sample calculation is desired, then no value for \code{y} should be specified.}
  \item{data}{If \code{x} is a formula, then \code{data} is an optional argument specifying data frame containing the variables in the formula.}
  \item{method}{Which version of the d statistic should we calculate? Possible values are \code{"pooled"} (the default), \code{"x.sd"}, \code{"y.sd"}, \code{"corrected"}, \code{"raw"}, \code{"paired"} and \code{"unequal"}. See below for specifics. }
  \item{mu}{The "null" value against which the effect size should be measured. This is almost always 0 (the default), so this argument is rarely specified.}
}

\details{The \code{cohensD} function calculates the Cohen's d measure of effect size in one of several different formats. The function is intended to be called in one of two different ways, mirroring the \code{t.test} function. That is, the first input argument \code{x} is a formula, then a command of the form \code{cohensD(x = outcome~group, data = data.frame)} is expected, whereas if \code{x} is a numeric variable, then a command of the form \code{cohensD(x = group1, y = group2)} is expected. Note that \code{cohensD} is not a generic function.

The \code{method} argument allows the user to select one of several different variants of Cohen's d. Assuming that the original t-test for which an effect size is desired was an independent samples t-test (i.e., not one sample or paired samples t-test), then there are several possibilities for how the normalising term (i.e., the standard deviation estimate) in Cohen's d should be calculated. The most commonly used method is to use the same pooled standard deviation estimate that is used in a Student t-test (\code{method = "pooled"}, the default). If \code{method = "raw"} is used, then the same pooled standard deviation estimate is used, except that the sample standard deviation is used (divide by N) rather than the unbiased estimate of the population standard deviation (divide by N-2). Alternatively, there may be reasons to use only one of the two groups to estimate the standard deviation. To do so, use \code{method = "x.sd"} to select the \code{x} variable, or the first group listed in the grouping factor; and \code{method = "y.sd"} to normalise by \code{y}, or the second group listed in the grouping factor. The last of the "Student t-test" based measures is the unbiased estimator of d (\code{method = "corrected"}), which multiplies the "pooled" version by (N-3)/(N-2.25). 

For other versions of the t-test, there are two possibilities implemented. If the original t-test did not make a homogeneity of variance assumption, as per the Welch test, the normalising term should mirror the Welch test (\code{method = "unequal"}). Or, if the original t-test was a paired samples t-test, and the effect size desired is intended to be based on the standard deviation of the differences, then \code{method = "paired"} should be used. 

The last argument to \code{cohensD} is \code{mu}, which represents an offset parameter. Specifically, it gives the value such that \code{mean(x) - mean(y) = mu} would be expected to correspond to an effect size estimate of 0. In most cases we would expect that \code{mu = 0} (the default), but not always.
}
\value{Numeric variable containing the effect size, d. Note that it does not show the direction of the effect, only the magnitude. That is, the value of d returned by the function is always positive or zero.
}


\references{
Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed.). Hillsdale, NJ: Lawrence Erlbaum Associates.
}
\author{
Daniel Navarro
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{t.test}}
}
\examples{
# calculate Cohen's d for two independent samples:
gradesA <- c(55, 65, 65, 68, 70) # 5 students with teacher A
gradesB <- c(56, 60, 62, 66)     # 4 students with teacher B
cohensD(gradesA, gradesB)

# calculate Cohen's d for the same data, described differently:
grade <- c(55, 65, 65, 68, 70, 56, 60, 62, 66) # grades for all students
teacher <- c("A", "A", "A", "A", "A", "B", "B", "B", "B") # teacher for each student
cohensD(grade ~ teacher)

# calculate Cohen's d for two paired samples:
pre  <- c(100, 122, 97, 25, 274) # a pre-treatment measure for 5 cases
post <- c(104, 125, 99, 29, 277) # the post-treatment measure for the same 5 cases
cohensD(pre, post, method = "paired") # ... explicitly indicate that it's paired, or else
cohensD(post - pre)  # ... do a "single-sample" calculation on the difference

# support for data frames:
exams <- data.frame(grade, teacher)
cohensD(exams$grade ~ exams$teacher)    # using $ 
cohensD(grade ~ teacher, data = exams)  # using the 'data' argument 

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
