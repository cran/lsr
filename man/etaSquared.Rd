\name{etaSquared}

\alias{etaSquared}

\title{Effect size calculations for ANOVAs}

\description{
Calculates eta-squared and partial eta-squared
}

\usage{
	etaSquared( x )
}

\arguments{
  \item{x}{An analysis of variance (aov) object.}
}

\details{
Calculates the eta-squared and partial eta-squared measures of effect size that are commonly used in analysis of variance. The input to the function should be the analysis of variance object itself.
}

\value{
The output is an M x 2 matrix. Each of the M rows corresponds to one of the terms in the ANOVA (e.g., main effect 1, main effect 2, interaction, etc), and each of the columns corresponds to a different measure of effect size. Column 1 contains the eta-squared values, and column 2 contains partial eta-squared values. 
}

\section{Warning}{This is an early version of the function, included only because of a looming teaching deadline, and has not been tested fully beyond the cases that I need to teach. Future versions may have somewhat different behaviour: e.g., extension to handle the output of multistratum analysis of variance objects.}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{Daniel Navarro}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}

\examples{
#### Example 1: one-way ANOVA ####

outcome <- c( 1.4,2.1,3.0,2.1,3.2,4.7,3.5,4.5,5.4 )  # data
treatment1 <- factor( c( 1,1,1,2,2,2,3,3,3 ))        # grouping variable
anova1 <- aov( outcome ~ treatment1 )                # run the ANOVA
summary( anova1 )                                    # print the ANOVA table
#
#             Df Sum Sq Mean Sq F value  Pr(>F)  
# treatment1   2 7.9356  3.9678  3.6626 0.09129 .
# Residuals    6 6.5000  1.0833
#
etaSquared( anova1 )                                 # effect size                            
#
#             eta.sq partial.eta.sq
# treatment1     0.5            0.5
#

#### Example 2: two-way ANOVA ####

treatment2 <- factor( c( 1,2,3,1,2,3,1,2,3 ))      # second grouping variable
anova2 <- aov( outcome ~ treatment1 + treatment2 ) # run the ANOVA
summary( anova2 )                                  # print the ANOVA table
#
#             Df Sum Sq Mean Sq F value   Pr(>F)   
# treatment1   2 7.9356  3.9678  55.797 0.001197 **
# treatment2   2 6.2156  3.1078  43.703 0.001915 **
# Residuals    4 0.2844  0.0711
#
etaSquared( anova2 )                               # effect size
#
#                eta.sq partial.eta.sq
# treatment1  0.5497229      0.9653961
# treatment2  0.4305727      0.9562393
#

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
