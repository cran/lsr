

#' Chi-square test against specified probabilities
#'
#' @description Convenience function that runs a chi-square goodness of fit
#' test against specified probabilities. This is a wrapper function intended
#' to be used for pedagogical purposes only.
#'
#' @param x Factor variable containing the raw outcomes.
#' @param p Numeric variable containing the null-hypothesis probabilities (default = all outcomes equally likely)
#'
#' @details The \code{goodnessOfFitTest} function runs the chi-square
#' goodness of fit test of the hypothesis that the outcomes in the factor
#' \code{x} were generated according to the probabilities in the vector
#' \code{p}. The probability vector \code{p} must be a numeric variable
#' of length \code{nlevels(x)}. If no probabilities are specified, all
#' outcomes are assumed to be equally likely.
#'
#' @return An object of class 'gofTest'. When printed, the output is
#' organised into four short sections. The first section lists the name
#' of the test and the variables included. The second lists the null and
#' alternative hypotheses for the test. The third shows the observed
#' frequency table, the expected frequency table under the null hypothesis,
#' and the probabilities specified by the null. The fourth prints out the
#' test results.
#'
#' @export
#'
#' @seealso
#' \code{\link{chisq.test}},
#' \code{\link{associationTest}},
#' \code{\link{cramersV}}
#'
#' @examples
#' # raw data
#' gender <- factor(
#'   c( "male","male","male","male","female","female",
#'      "female","male","male","male" ))
#'
#' # goodness of fit test against the hypothesis that males and
#' # females occur with equal frequency
#' goodnessOfFitTest( gender )
#'
#' # goodness of fit test against the hypothesis that males appear
#' # with probability .6 and females with probability .4.
#' goodnessOfFitTest( gender, p=c(.4,.6) )
#' goodnessOfFitTest( gender, p=c(female=.4,male=.6) )
#' goodnessOfFitTest( gender, p=c(male=.6,female=.4) )
#'
goodnessOfFitTest <- function( x, p=NULL ) {

  # check if x is missing
  if( missing(x) ) { stop( '"x" argument is missing, with no default')}

  # make sure x is a factor
  if( !methods::is(x,"factor") ) {
    stop( "input argument 'x' must be a factor" )
  }

  # check for missing data & print warning if needed
  missingData <- is.na( x )
  x <- x[!missingData]
  if(sum(missingData) > 0) {
    warning( paste(sum(missingData)), " case(s) removed due to missingness" )
  }

  # make sure x has at least two levels
  if( nlevels(x) <= 1 ) {
    stop( "factor 'x' must be have at least two levels" )
  }

  # either set the default probabilities...
  if( is.null(p) ) {
    p <- rep.int( 1, nlevels(x) ) / nlevels(x)
    names(p) <- levels( x )
  } else { # or check that the user supplied probabilites make sense...
    if( !methods::is(p,"numeric") | any(is.na(p)) ) stop( "'p' must be a numeric vector of probabilities" )
    if( length(p) != nlevels(x) ) stop( "'p' contains the wrong number of elements" )
    if( abs(sum(p) - 1)>10^-10 ) {
      warning( "probabilities in 'p' do not add up to 1. rescaled values used")
      p <- p / sum(p)
    }
    if( !is.null( names(p)) ) { # if p has names, better check
      if( all ( sort(names(p)) == sort(levels(x)) )) { # if they match
        p <- p[levels(x)] # sort them to match level ordering
      } else {
        warning( "'p' has names that differ from the levels of 'x'")
      }
    }
  }

  # tabulate x
  f <- table(x)

  # run the corresponding chi-square test: suppressing the warning,
  # replacing it with our own if it exists
  old.warn <- options(warn=2) # convert warnings to errors
  htest <- try( stats::chisq.test( x=f, p=p ), silent=TRUE ) # try the test
  need.warning <- FALSE # assume warning

  # check for failures:
  if( class(htest) == "try-error" ) {

    # handle the case when the "error" was the warning we're catching
    need.warning <- length( grep("Chi-squared approximation may be incorrect",htest )) > 0
    if( need.warning ) { # if it was a ties problem...
      options(warn=-1) # suppress warnings completely for the next run...
    } else { # if not...
      options( old.warn ) # reset warning state and let R throw what it likes...
    }

    # now run the test again & compute effect size
    htest <- stats::chisq.test( x=f, p=p )
  }
  options( old.warn )   # reset warnings

  # get the variable name if it exists
  outcome <- match.call()[2] # get the x argument from the call
  outcome <- as.character( outcome )

  # create output structure
  out <- list(
    outcome = outcome,
    p = p,
    observed = htest$observed,
    expected = htest$expected,
    difference = htest$observed - htest$expected,
    statistic = htest$statistic,
    df = htest$parameter,
    p.value = htest$p.value,
    warn = need.warning
  )
  class( out ) <- "gofTest"

  # throw the warning if needed
  if( need.warning ) {
    warning( "Expected frequencies too small: chi-squared approximation may be incorrect")
  }

  return(out)

}


#' Print method for lsr goodness-of-fit tests
#'
#' @param x An object of class 'gofTest'
#' @param ... For consistency with the generic (unused)
#'
#' @return Invisibly returns the original object
#' @export
print.gofTest <- function(x, ...) {

  # print the name of the test
  cat("\n     Chi-square test against specified probabilities\n\n")

  # print the data variable
  cat( "Data variable:  ", x$outcome, "\n\n" )

  # print the hypotheses being tested
  cat( "Hypotheses: \n")
  cat( "   null:        true probabilities are as specified\n")
  cat( "   alternative: true probabilities differ from those specified\n")
  cat("\n")

  # print observed against expected
  descriptives <- cbind( x$observed, x$expected, x$p )
  colnames(descriptives) <- c("observed freq.",  "expected freq.", "specified prob." )
  cat( "Descriptives: \n")
  print(descriptives)
  cat( "\n")

  # print test statistics
  nDigits <- 3
  cat( "Test results: \n")
  cat( "   X-squared statistic: ", round( x$statistic, nDigits ), "\n" )
  cat( "   degrees of freedom: ", round( x$df, nDigits ), "\n" )
  pp <- ifelse( x$p.value < .001, "<.001", round( x$p.value, nDigits ))
  cat( "   p-value: ", pp, "\n")
  if( x$warn) cat( "   warning: expected frequencies too small, results may be inaccurate\n")
  cat( "\n")


  return( invisible(x) )

}

