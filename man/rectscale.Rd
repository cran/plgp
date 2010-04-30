\name{rectscale}
\alias{rectscale}
\alias{rectunscale}

\title{
  Un/Scale data in a bounding rectangle
}
\description{
  Scale data lying in an arbitrary rectangle to lie in
  the unit rectangle, and back again
}
\usage{
rectscale(X, rect)
rectunscale(X, rect)
}

\arguments{
  \item{X}{
    a \code{matrix} or \code{data.frame} of real-valued covariates
  }
  \item{rect}{
    a \code{matrix} describing a bounding rectangle for \code{X}
    with 2 columns and \code{ncol(X)} rows
  }
}
\value{
  a \code{matrix} or \code{data.frame} with the same dimensions as
  \code{X} scaled or un-scaled as appropriate
}
\author{
  Robert B. Gramacy, \email{bobby@statslab.cam.ac.uk}
}
\references{
  \url{http://www.statslab.cam.ac.uk/~bobby/plgp.html}
}
\examples{
X <- matrix(runif(10, 1, 3), ncol=2)
rect <- rbind(c(1,3), c(1,3))
Xs <- rectscale(X, rect)
rectunscale(Xs, rect)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ utilities }
