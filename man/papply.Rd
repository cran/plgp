\name{papply}
\alias{papply}
\title{
  Extending apply to particles
}
\description{
  Applies a user-specified function to each particle contained in the
  global variables \code{peach} and \code{pall}, collecting the
  output in a \code{\link{data.frame}}
}
\usage{
papply(fun, verb = 1, pre = "", ...)
}

\arguments{
  \item{fun}{
    a user-defined function which which takes a particle as its first
    input; the output of \code{fun} should be a vector, \code{matrix}
    or \code{data.frame}
  }
  \item{verb}{
    a scalar \code{logical} indicating whether progress statements
    should be printed to the screen
  }
  \item{pre}{
    an optional \code{character} prefix used in the progress print
    statements; ignored if \code{verb = 0}
  }
  \item{\dots}{
    these ellipses arguments are used to pass extra optional
    arguments to the user-supplied function \code{fun}
  }
}
\details{
  This is a extension to the built-in \code{\link{apply}} family of
  function to particles, intended to be used with the particles created
  by \code{\link{PL}}.  Perhaps the most common use of this function is
  in obtaining samples form the posterior predictive distribution, i.e.,
  with the user supplied \code{fun = \link{pred.GP}}

  The particles applied over must be present in the global variables
  \code{pall}, containing sufficient information common to all
  particles, \code{peach}, containing sufficient information
  particular to each particle, as constructed by \code{\link{PL}}
}
\value{
  Returns a data frame with the collected output of the user-specified
  function \code{fun}
}
\references{
  Carvalho, C., Johannes, M., Lopes, H., and Polson, N. (2008).
  \dQuote{Particle Learning and Smoothing.}
  Discussion Paper 2008-32, Duke University Dept. of Statistical
  Science.

  \url{http://faculty.chicagobooth.edu/robert.gramacy/plgp.html}
}
\author{
Robert B. Gramacy, \email{rbgramacy@chicagobooth.edu}
}
\seealso{
  \code{\link{PL}}, \code{\link{pred.GP}}
}
\examples{
## See the demos via demo(package="plgp") and the examples
## section of ?plgp
}
\keyword{iteration}
\keyword{methods}
