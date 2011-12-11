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

This function is currently just a wrapper function for \code{\link{pairwise.t.test}}, and it only works for one-way ANOVA, but this will change in future versions. Additionally, the structure of the output will change so that the output is closer in style to \code{TukeyHSD}.
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
#
#		  outcome group
#		1       1     a
#		2       2     a
#		3       3     a
#		4       2     b
#		5       3     b
#		6       4     b
#		7       5     c
#		8       6     c
#		9       7     c
#

# run the ANOVA and print out the ANOVA table:
anova1 <- aov( outcome ~ group, data = dataset )
summary(anova1)
#
#             Df Sum Sq Mean Sq F value   Pr(>F)   
# group        2     26      13      13 0.006592 **
# Residuals    6      6       1
#

# The F-test is significant, and we had no theoretical reason
# think that any particular group difference would be the  
# source of the effect. So we run post hoc t-tests comparing 
# all pairs of groups, making adjustments to control the 
# overall Type I error rate. Since Bonferroni (and Holm) 
# methods are among the easiest to understand, and are typically
# the first methods taught in undergraduate classes, we'll use
# those. 
#
# Currently, the following two commands are equivalent:

posthocPairwiseT( anova1 )
pairwise.t.test( dataset$outcome, dataset$group )
#
# 	 Pairwise comparisons using t tests with pooled SD 
#
# data:  dataset$outcome and dataset$group 
#
#   a      b     
# b 0.2666 -     
# c 0.0081 0.0208
# 
# P value adjustment method: holm
# 

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
