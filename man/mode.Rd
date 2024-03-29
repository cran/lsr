% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maxFreq.R, R/modeOf.R
\name{maxFreq}
\alias{maxFreq}
\alias{mode}
\alias{modeOf}
\title{Sample mode}
\usage{
maxFreq(x, na.rm = TRUE)

modeOf(x, na.rm = TRUE)
}
\arguments{
\item{x}{A vector containing the observations.}

\item{na.rm}{Logical value indicating whether NA values should be removed.}
}
\value{
The \code{modeOf} function returns the mode of \code{x}. If there
are ties, it returns a vector containing all values of \code{x} that have
the modal frequency. The \code{maxFreq} function returns the modal
frequency as a numeric value.
}
\description{
Calculate the mode of a sample: both modal value(s) and
the corresponding frequency
}
\details{
These two functions can be used to calculate the mode (most
frequently observed value) of a sample, and the actual frequency of the
modal value. The only complication is in respect to missing data. If
\code{na.rm = FALSE}, then there are multiple possibilities for how to
calculate the mode. One possibility is to treat \code{NA} as another
possible value for the elements of \code{x}, and therefore if \code{NA}
is more frequent than any other value, then \code{NA} is the mode; and
the modal frequency is equal to the number of missing values. This is
the version that is currently implemented.

Another possibility is to treat \code{NA} as meaning "true value unknown",
and to the mode of \code{x} is itself known only if the number of missing
values is small enough that -- regardless of what value they have -- they
cannot alter the sample mode. For instance, if \code{x} were
\code{c(1,1,1,1,2,2,NA)}, we know that the mode of \code{x} is \code{1}
regardless of what the true value is for the one missing datum; and we
know that the modal frequency is between 4 and 5. This is also a valid
interpretation, depending on what precisely it is the user wants, but
is not currently implemented.

Because of the ambiguity of how \code{na.rm = FALSE} should be interpreted,
the default value has been set to \code{na.rm = TRUE}, which differs from
the default value used elsewhere in the package.
}
\examples{
# simple example
eyes <- c("green","green","brown","brown","blue")
modeOf(eyes)
maxFreq(eyes)

# vector with missing data
eyes <- c("green","green","brown","brown","blue",NA,NA,NA)

# returns NA as the modal value.
modeOf(eyes, na.rm = FALSE)
maxFreq(eyes, na.rm = FALSE)

# returns c("green", "brown") as the modes, as before
modeOf(eyes, na.rm = TRUE)
maxFreq(eyes, na.rm = TRUE)
}
\seealso{
\code{\link{mean}},
\code{\link{median}},
\code{\link{table}}
}
