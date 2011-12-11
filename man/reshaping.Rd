\name{reshaping}
\alias{wideToLong}
\alias{longToWide}
\alias{wideToMV}
\alias{wideRM}
\alias{longRM}


\title{
Reshaping functions for repeated measures data
}

\description{
Functions to switch between wide, long, and multivariate form of a data frame for a repeated measures study design 
}

\usage{
# reshaping functions:
wideToLong( data, rms = NULL, ... )
longToWide( data, rms = NULL, ... )
wideToMV( data, rms = NULL, ... )

# functions to create a repeated measures structure:
wideRM( data, treatments = NULL, sep = "_" )
longRM( data, treatments, measures, between, sep = "_" )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{The data frame.}
  \item{rms}{A repeated measures structure (see below).}
  \item{...}{Additional arguments to be passed to the repeated measures structure functions.}
  \item{treatments}{Names of the within-subjects treatment variables.}
  \item{measures}{Names of the measured (outcome) variables.}
  \item{between}{Names of the between-subjects variables.}
  \item{sep}{A separator character (or string).}
}

\details{
These functions serve two  purposes: firstly, the \code{wideRM} and \code{longRM} functions can be used to extract the repeated measures structure (see below) from a data frame that contains repeated measurements. Secondly, once a repeated measures structure has been specified, the other three functions (\code{wideToLong}, \code{longToWide} and \code{wideToMV}) can be used to easily switch between three standard formats for a data set containing repeated measurements.

All three of the reshaping functions require a repeated measures structure for the optional argument \code{rms}. If a correctly specified repeated measures structure exists, then the reshaping functions impose no other constraints: the commands \code{wideToLong(data,rms)}, \code{longToWide(data,rms)} and \code{wideToMV(data,rms)} all work without additional arguments or any other requirements. If the \code{rms} structure is not specified, then the reshaping function will attempt to construct one: the additional arguments in \code{...} are passed to the \code{wideRM} function (from \code{wideToLong} and \code{wideToMV}) or to the \code{longRM} (from \code{longToWide}). 

\strong{Long form and wide form}

In the wide format each row corresponds to single subject (or other experimental unit) and each column is a variable. If, as is typical of repeated measures designs, each participant is measured under multiple conditions using the same outcome variable, then each measurement would constitute a different variable. For instance, if the experiment used 10 subjects, and measured variable \code{outcome} in two conditions \code{cond1} and \code{cond2}, then the wide form data frame would have 10 rows, and would have separate columns for each condition, likely named \code{outcome_cond1} and \code{outcome_cond2}. One would typically also have an \code{id} variable that uniquely identifies each subject:

%\preformatted{
%    id outcome_cond1 outcome_cond2
%  1  1           415           455
%  2  2           500           532
%  3  3           478           499
%  4  4           550           602 
%}

While a wide form data frame works on a one row per subject basis, a long form data frame has one row per subject-condition combination. Instead of having separate \code{outcome} variables for each condition, there is only a single \code{outcome} variable, plus a separate (factor) variable indicating the condition in which each outcome was observed. Note that if there are multiple treatments (e.g., \code{condition} and \code{block}), then the long form of the data frame would contain two within subject factors. Similarly, if there are multiple response variables, then there can be two outcomes, e.g., \code{outcome1} and \code{outcome2}.

%The equivalent long form data frame to the wide form data set shown above is:
%\preformatted{
%    id outcome condition 
%  1  1     415     cond1 
%  2  1     455     cond2
%  3  2     500     cond1
%  4  2     532     cond2
%  5  3     478     cond1
%  6  3     499     cond2
%  7  4     550     cond1
%  8  4     602     cond2
%}



The multivariate form, which is less frequently needed in an introductory statistics setting and is not in fact a data frame, is described at the end of the section. 

\strong{Repeated measures structures}

A repeated measures structure, as used in these functions, is a list with two components. If the repeated measures structure is \code{rms}, then \code{rms$between} is a character vector listing the names of the between subjects variables. In contrast, \code{rms$within} is a data frame that describes both the long form and wide form names for the repeated measurements. Specifically, \code{rms$within} has one row per repeated-measures variable (in the wide form data frame). The first column of \code{rms$within} is always called \code{wide.name}, and contains the name of the variable as it appears in the wide form of the data set. The second column is called \code{measure} and contains the name of the outcome variable being measured. Each subsequent column represents one of the repeated-measures variable: the name of the variable is the name of the repeated-measures factor, and the values of the variable indicate which condition is represented by the wide form variable. See the Examples section for more detail. 


%To make this more concrete, the value of \code{rms$between$} for the example discussed above would just be a single character string (i.e. \code{"id"}) since there is only one between subjects variable, and the value of \code{rms$within} would be the following data frame:
%\verb{ 
%        wide.name measure condition
%  1 outcome_cond1 outcome     cond1
%  2 outcome_cond2 outcome     cond2
%}

A repeated measures structure can be specified manually if need be, but in general it is easier to use the \code{wideRM} function to extract it from a wide form data frame, or to use the \code{longRM} function to extract it from the long form data. 

The \code{wideRM} function is simple to use, but quite restrictive. It relies entirely on the assumption that the variable names in the wide form data frame follow a strict naming convention. Specifically, for a data set with two within subjects factors, the within subjects variables need to be named according to the following scheme:

\code{[measurement][sep][treatment1][sep][treatment2]} 

where \code{sep} is a separator string that does not appear anywhere else in the variable names (including the between subjects variable names). By default, it is assumed that \code{sep = "_"}. For example, suppose we have an experiment with two within-subjects factors (e.g., \code{block} and \code{difficulty}), and two outcome measures (e.g., mean response time \code{MRT} and proportion of correct responses \code{PC}). Then the wide form data frame would contain 8 within-subject variables, and a valid naming scheme for these would be to name them as follows:

\code{ 
  MRT_block1_easy,
  MRT_block1_hard,
  MRT_block2_easy,
  MRT_block2_hard,
  PC_block1_easy,
  PC_block1_hard,
  PC_block2_easy,
  PC_block2_hard
}

As long as the within-subjects variables follow this kind of format, and that the separator string (which in this case is the default value of \code{sep = "_"}) does not appear in the names of the between-subjects variables, then the variable names uniquely specify everything about the repeated measures structure \code{rms} except for the names of within-subjects factors. To specify these names, use the \code{treatments} argument. The names should appear in the same order that the corresponding values appear in the wide form variable names. In this case, a value of \code{treatment = c("block","difficulty")} would be appropriate. If no value for \code{treatment} is given, then the within-subjects factors are labelled \code{treatment.1}, \code{treatment.2} and so on.

Constructing a repeated measures structure from the long form of a data frame is also possible, using the \code{longRM} function. In addition to (optionally) specifying the separator string \code{sep}, in this case the user needs to specify which of the variables in the long form data frame are the \code{between} subjects variables, which ones are the \code{treatments} and which ones are the \code{measures}. For instance, to continue the example used when discussing \code{wideRM} above, the appropriate values might be something like \code{between = "id"}, \code{treatments = c("block","difficulty")} and \code{measures = c("MRT","PC")}. 

\strong{The multivariate form}

The final format that is supported by these functions is the multivariate form. This is essentially equivalent to the wide form of the data frame, except for the fact that all wide form variables that correspond to the same measured variable (e.g. all of the \code{MRT} variables) are grouped together into a matrix. As a consequence, the multivariate form is a list rather than a data frame: each between subjects variable is a vector element to the list, and each within subjects variable is a matrix. This format is used much less frequently by beginners, and so is not documented in as much detail. See the Examples section for more information.

\strong{Final remarks}

Overall, the reshaping framework is less flexible than the "\code{melt} and \code{cast}" approach provided by the \pkg{reshape} package; it is similar to but less flexible than the \code{\link{reshape}} function. This is deliberate: by restricting the range of possibilities considerably it allows for a much simpler syntax, and has the desirable side-effect that there is a simple relationship between the repeated measure structure produced by \code{wideRM} and \code{longRM} and the \code{idesign} matrix that is required to run repeated measures analysis of variance models using the \code{Anova} functions in the \pkg{car} package.
}




\value{The reshaping functions \code{wideToLong} and \code{longToWide} return data frames in the appropriate format, whereas the \code{wideToMV} returns a list. The other two functions \code{wideRM} and \code{longRM} both return repeated measures structures, lists containing two elements. 
}

%\references{
%% ~put references to the literature/web site here ~
%}

\author{
Daniel Navarro
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{reshape}}, \code{\link{stack}}, \code{cast} and \code{melt} (in the \pkg{reshape} package)
}


\examples{
##### Example 1 : single repeated measure design #####

# Outcome measure is mean response time (MRT), measured in two conditions
# with 4 participants. All participants participate in both conditions. 

MRT_cond1 <- c( 415,500,478,550 )   # response times in condition 1
MRT_cond2 <- c( 455,532,499,602 )   # response times in condition 2
id <- 1:4                           # id variable

# Note that all repeated measures variable names are formatted in the
# "right" format: the name of the outcome variable (MRT), followed by a 
# separator character (_) and then the value on the repeated measures
# factor (e.g. cond1).  Similarly, note that the separator character
# does not appear anywhere else in the names; including the names of
# the between subject variable (i.e. id). 

expt <- data.frame( id, MRT_cond1, MRT_cond2 ) # convert to data frame

# This is what the expt data frame looks like:
#
#   id MRT_cond1 MRT_cond2
# 1  1       415       455
# 2  2       500       532
# 3  3       478       499
# 4  4       550       602
#
# This is a standard "wide" form data frame. 

# --- Example 1.1 --- Extracting the repeated measures structure from  
# the wide form data frame is straightforward in this case, because all
# the variable names match the defaults

wideRM( expt )

# This is what the output looks like:
# 
# $within
#   wide.name measure treatment.1
# 1 MRT_cond1     MRT       cond1
# 2 MRT_cond2     MRT       cond2
#
# $between
# [1] "id"
#

# The equivalent command with all input arguments specified manaully 
# looks like this:

wideRM( data = expt,                  # the data frame
        treatments = "treatment.1",   # name of the within-subjects treatment
        sep = "_"                     # the separator string used in the names
)


# --- Example 1.2 --- Converting the wide form data frame to the long form 
# data frame is equally simple, again because the variable names have been
# specified in the "right" format:

wideToLong( expt )

# Here's the output:
#
#   id MRT treatment.1
# 1  1 415       cond1
# 2  1 455       cond2
# 3  2 500       cond1
# 4  2 532       cond2
# 5  3 478       cond1
# 6  3 499       cond2
# 7  4 550       cond1
# 8  4 602       cond2

# And the equivalent command with all arguments specified manually:

wideToLong( data = expt,          # the data frame
            rms = wideRM( expt )  # the repeated measures structure 
)

# --- Example 1.3 --- Conversion from long form to wide form is 
# straightforward as long as we have the repeated measures structure
# available to us. However, constructing the repeated measures
# structure from a long form data frame is slightly more tedious. 
# First, we need a long form data frame:

expt2 <- wideToLong( expt )

# This is the same long form data frame above. We can create the 
# repeated measures structure as follows:

expt2.rms <- longRM( data = expt2,                # data frame
                     treatments = "treatment.1",  # within subjects treatments
                     measures = "MRT",            # measured variables (outcomes)
                     between = "id",              # between subjects variables
                     sep = "_"                    # separator character
                   )

# This produces the exact same repeated measures structure that we
# obtained from the wideRM(expt) command earlier. 

# Now that we have a repeated measures structure specified, the 
# reshaping is straightforward:

longToWide( data = expt2, rms = expt2.rms )

# The output here is identical to the wide form data frame that we 
# started with.


# --- Example 1.4 --- Conversion from wide form to a "multivariate" form.
# This is useful for multivariate linear models. As before, it's easy in 
# this case because the names in the wide form data are structured the 
# way we need it:

wideToMV( expt )
wideToMV( data = expt, rm = wideRM(expt) )  # equivalent command

# Here's the output: 
# 
# $MRT
#      MRT_cond1 MRT_cond2
# [1,]       415       455
# [2,]       500       532
# [3,]       478       499
# [4,]       550       602
# 
# $id
# [1] 1 2 3 4
# 


##### Example 2 : two treatments and two measures #####

# A more complex, but more realistic, version of the experiment. Again, we have only
# four participants, but now we have two different outcome measures, mean response
# time (MRT) and the proportion of correct responses (PC). Additionally, we have two
# different repeated measures variables. As before, we have the experimental condition
# (cond1, cond2), but this time each participant does both conditions on two different
# days (day1, day2). Finally, we have multiple between-subject variables too, namely
# id and gender.

# response times across both conditions and both days:
MRT_cond1_day1 <- c( 415,500,478,550 )
MRT_cond2_day1 <- c( 455,532,499,602 )
MRT_cond1_day2 <- c( 400,490,468,502 )
MRT_cond2_day2 <- c( 450,518,474,588 )

# proportion of correct reponses in both conditions and days:
PC_cond1_day1 <- c( 79,83,91,75 )
PC_cond2_day1 <- c( 82,86,90,78 )
PC_cond1_day2 <- c( 88,92,98,89 )
PC_cond2_day2 <- c( 93,97,100,95 )

# between subjects variables
id <- 1:4
gender <- factor( c("male","male","female","female") )

# create wide form data frame
expt3 <- data.frame(  id, gender, 
                      MRT_cond1_day1, MRT_cond1_day2, MRT_cond2_day1, MRT_cond2_day2, 
                      PC_cond1_day1, PC_cond1_day2, PC_cond2_day1, PC_cond2_day2 
                    )

# Here's the wide form data frame:
#
#   id gender MRT_cond1_day1 MRT_cond1_day2 MRT_cond2_day1 MRT_cond2_day2 PC_cond1_day1 PC_cond1_day2 PC_cond2_day1 PC_cond2_day2
# 1  1   male            415            400            455            450            79            88            82            93
# 2  2   male            500            490            532            518            83            92            86            97
# 3  3 female            478            468            499            474            91            98            90           100
# 4  4 female            550            502            602            588            75            89            78            95


# Extracting the repeated measures structure from the variable names:
wideRM( expt3 )

# Output:
# 
# $within
#        wide.name measure treatment.1 treatment.2
# 1 MRT_cond1_day1     MRT       cond1        day1
# 2 MRT_cond1_day2     MRT       cond1        day2
# 3 MRT_cond2_day1     MRT       cond2        day1
# 4 MRT_cond2_day2     MRT       cond2        day2
# 5  PC_cond1_day1      PC       cond1        day1
# 6  PC_cond1_day2      PC       cond1        day2
# 7  PC_cond2_day1      PC       cond2        day1
# 8  PC_cond2_day2      PC       cond2        day2
#
# $between
# [1] "id"     "gender"
#

# Conversion to long form:
wideToLong( expt3 )

# Output: 
#
#   id gender MRT  PC treatment.1 treatment.2
# 1   1   male 415  79       cond1        day1
# 2   1   male 400  88       cond1        day2
# 3   1   male 455  82       cond2        day1
# 4   1   male 450  93       cond2        day2
# 5   2   male 500  83       cond1        day1
# 6   2   male 490  92       cond1        day2
# 7   2   male 532  86       cond2        day1
# 8   2   male 518  97       cond2        day2
# 9   3 female 478  91       cond1        day1
# 10  3 female 468  98       cond1        day2
# 11  3 female 499  90       cond2        day1
# 12  3 female 474 100       cond2        day2
# 13  4 female 550  75       cond1        day1
# 14  4 female 502  89       cond1        day2
# 15  4 female 602  78       cond2        day1
# 16  4 female 588  95       cond2        day2
#

# Conversion to multivariate form:
wideToMV( expt3 )

# Output:
#
# $MRT
#      MRT_cond1_day1 MRT_cond1_day2 MRT_cond2_day1 MRT_cond2_day2
# [1,]            415            400            455            450
# [2,]            500            490            532            518
# [3,]            478            468            499            474
# [4,]            550            502            602            588
#
# $PC
#      PC_cond1_day1 PC_cond1_day2 PC_cond2_day1 PC_cond2_day2
# [1,]            79            88            82            93
# [2,]            83            92            86            97
# [3,]            91            98            90           100
# [4,]            75            89            78            95
#
# $id
# [1] 1 2 3 4
#
# $gender
# [1] male   male   female female
# Levels: female male
#

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
