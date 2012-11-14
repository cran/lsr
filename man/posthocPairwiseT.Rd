\name{posthocPairwiseT}
\alias{posthocPairwiseT}

\title{
Post-hoc pairwise t-tests for ANOVA
}

\description{
Performs pairwise t-tests for an analysis of variance, making corrections for multiple comparisons. 
}

\usage{
posthocPairwiseT(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An \code{aov} objects}
  \item{...}{Arguments to be passed to \code{pairwise.t.test}}
}

\details{
The intention behind this function is to allow users to use simple tools for multiple corrections (e.g., Bonferroni, Holm) as post hoc corrections in an ANOVA context, using the fitted model object (i.e., an \code{aov} object) as the input. The reason for including this function is that Tukey / Scheffe methods for constructing simultaneous confidence intervals (as per \code{\link{TukeyHSD}}) are not often discussed in the context of an introductory class, and the more powerful tools provided by the \code{multcomp} package are not appropriate for students just beginning to learn statistics.

This function is currently just a wrapper function for \code{\link{pairwise.t.test}}, and it only works for one-way ANOVA, but this may change in future versions. 
}

\value{As per \code{pairwise.t.test}}

\section{Warning}{This function is under development. The internals to the function will almost certainly change: in particular the \code{p.adjust.method} argument from the \code{pairwise.t.test} will very likely change to \code{method}, matching the \code{p.adjust} function instead.}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Daniel Navarro
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pairwise.t.test}}, \code{\link{TukeyHSD}}, \code{glht} (in the \code{multcomp} package)
}


\examples{# create the data set to analyse: 	
dataset <- data.frame( 
	        outcome = c( 1,2,3, 2,3,4, 5,6,7 ),
            group = factor(c( "a","a","a", "b","b","b","c","c","c")) 
           )

# run the ANOVA and print out the ANOVA table:
anova1 <- aov( outcome ~ group, data = dataset )
summary(anova1)

# Currently, the following two commands are equivalent:
posthocPairwiseT( anova1 )
pairwise.t.test( dataset$outcome, dataset$group )


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
