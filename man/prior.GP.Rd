\name{prior.GP}
\alias{prior.GP}
\alias{prior.CGP}
\alias{prior.ConstGP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate priors for GP models
}
\description{
  Generate priors for Gaussian process (GP)
  regression, classification, or combined unknown constraint
  models
}
\usage{
prior.GP(m, cov = c("isotropic", "separable", "sim"))
prior.CGP(m, cov = c("isotropic", "separable", "sim"))
prior.ConstGP(m, cov.GP = c("isotropic", "separable", "sim"),
              cov.CGP = cov.GP)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
    positive scalar integer specifying the dimensionality of the input
    space
  }
  \item{cov}{
    whether to use an \code{"isotropic"} or \code{"separable"} power
    exponential correlation function with power 2 -- nugget included;
    a single index model (\code{"sim"}) capability is provided as \dQuote{beta}
    functionality; applies to both regression and classification GPs
  }
  \item{cov.GP}{
    specifies the covariance for the real-valued response in the
    combined unknown constraint GP model
  }
  \item{cov.CGP}{
   specifies the covariance for the categorical response in the
   combined unknown constraint GP model
  }
}
\details{
  These function generate a default prior object in the correct format
  for use with the other PL routines, e.g.,
  \code{\link{init.GP}} and \code{\link{pred.GP}}.  The object returned
  may be modified as necessary.

  The \code{\link{prior.ConstGP}} is essentially the combination
  of \code{\link{prior.GP}} and \code{\link{prior.CGP}}
  for regression and classification GP models, respectively
}
\value{
  a valid prior object for the appropriate GP model;

  By making the output \code{$drate} and/or \code{$grate}
  values negative causes the corresponding lengthscale \code{d}
  parameter(s) and nugget \code{d} parameter to be fixed at the
  reciprocal of their absolute values, respectively.  This effectively
  turns off inference for these values, and allows one to study the GP
  predictive distribution as a function of fixed values.  When both
  are fixed it is sensible to use only one particle (\code{P=1}, as an
  argument to \code{\link{PL}}) 
}
\references{

 Gramacy, R. and Polson, N. (2011).
  \dQuote{Particle learning of Gaussian process models for
    sequential design and optimization.}
  Journal of Computational and Graphical Statistics, 20(1), 
  pp. 102-118; arXiv:0909.5262

  Gramacy, R. and Lee, H. (2010).
  \dQuote{Optimization under unknown constraints}.
  \emph{Bayesian Statistics 9}, J. M. Bernardo, M. J. Bayarri,
  J. O. Berger, A. P. Dawid, D. Heckerman, A. F. M. Smith and M. West
  (Eds.); Oxford University Press

  Gramacy, R. (2020).
  \dQuote{Surrogates: Gaussian Process Modeling, Design and Optimization for the Applied Sciences}.
  Chapman Hall/CRC; \url{https://bobby.gramacy.com/surrogates/}

  \url{https://bobby.gramacy.com/r_packages/plgp/}
}
\author{
Robert B. Gramacy, \email{rbg@vt.edu}
}
\seealso{
  \code{\link{PL}}, \code{\link{lpredprob.GP}},
  \code{\link{propagate.GP}}, \code{\link{init.GP}},
  \code{\link{pred.GP}}
}
\examples{
## See the demos via demo(package="plgp") and the examples
## section of ?plgp
}

\keyword{models}
\keyword{regression}
\keyword{classif}
\keyword{methods}

