\name{pairedSamplesTTest}
\alias{pairedSamplesTTest}


\title{
Paired samples t-test
}

\description{
Convenience function that runs a paired samples t-test. This is a wrapper function intended to be used for pedagogical purposes only.
}

\usage{
pairedSamplesTTest( 
	 formula, data=NULL, id=NULL,
	 one.sided=FALSE, conf.level=.95 
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{Formula specifying the outcome and the groups (required).}
  \item{data}{Optional data frame containing the variables.}
  \item{id}{The name of the id variable (must be a character string).}
  \item{one.sided}{One sided or two sided hypothesis test (default = \code{FALSE})}
  \item{conf.level}{The confidence level for the confidence interval (default = .95).}
}

\details{
The \code{pairedSamplesTTest} function runs a paired-sample t-test, and prints the results in a format that is easier for novices to handle than the output of \code{t.test}. All the actual calculations are done by the \code{t.test} and \code{cohensD} functions. The data must be formatted as a "long form" data frame, as is typical of other functions that implement repeated measures analyses. The \code{formula} argument must be a two-sided formula of the form \code{outcome ~ group}, and the \code{id} argument must be a string specifying the name of the id variable. 

As with the \code{t.test} function, the default test is two sided, corresponding to a default value of \code{one.sided = FALSE}. To specify a one sided test, the \code{one.sided} argument must specify the name of the factor level that is hypothesised (under the alternative) to have the larger mean. For instance, if the outcome at "time2" is expected to be higher than at "time1", then the corresponding one sided test is specified by \code{one.sided = "time2"}. 
}

\value{An object of class 'TTest'. When printed, the output is organised into five short sections. The first section lists the name of the test and the variables included. The second provides means and standard deviations. The third states explicitly what the null and alternative hypotheses were. The fourth contains the test results: t-statistic, degrees of freedom and p-value. The final section includes the relevant confidence interval and an estimate of the effect size (i.e., Cohen's d)}

%\references{
%% ~put references to the literature/web site here ~
%}

\section{Note}{Although not "properly" documented, \code{pairedSamplesTTest} also supports an "lme4"-like method for specifying the entire model in the \code{formula} argument, removing the need to specify the \code{id} argument. To do so, the formula should be of the form \code{outcome ~ group + (id)} or \code{outcome ~ group + (1|id)}.}

\section{Warning}{This package is under development, and has been released only due to teaching constraints. Until this notice disappears from the help files, you should assume that everything in the package is subject to change. Backwards compatibility is NOT guaranteed. Functions may be deleted in future versions and new syntax may be inconsistent with earlier versions. For the moment at least, this package should be treated with extreme caution.}

\author{
Daniel Navarro
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{t.test}}, \code{\link{oneSampleTTest}}, \code{\link{independentSamplesTTest}}, \code{\link{cohensD}}
}


\examples{

df <- data.frame(
       id = factor( x=c(1, 1, 2, 2, 3, 3, 4, 4), 
                    labels=c("alice","bob","chris","diana") ), 
       time = factor( x=c(1,2,1,2,1,2,1,2), 
                      labels=c("time1","time2")), 
       wm = c(3, 4, 6, 6, 9, 12,7,9) 
)
	
	
pairedSamplesTTest( formula= wm~time, data=df, id="id" )

# one sided test
pairedSamplesTTest( formula= wm~time, data=df, id="id", one.sided="time2" )

# missing data because of NA values
df$wm[1] <- NA
pairedSamplesTTest( formula= wm~time, data=df, id="id" )

# missing data because of missing cases from the data frame
df <- df[-1,] 
pairedSamplesTTest( formula= wm~time, data=df, id="id" )


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
