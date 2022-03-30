\name{addpall.GP}
\alias{addpall.GP}
\alias{addpall.CGP}
\alias{addpall.ConstGP}

\title{
  Add data to pall
}
\description{
  Add sufficient 
  data common to all particles to the global \code{pall}
  variable, a mnemonic for \dQuote{particles-all}, for
  Gaussian process (GP)
  regression, classification, or combined unknown constraint
  models
}
\usage{
addpall.GP(Z)
addpall.CGP(Z)
addpall.ConstGP(Z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
    new observation(s) (usually the next one in \dQuote{time}) to add to
    the \code{pall} global variable 
  }
}
\details{
  All three functions add new \code{Z$x} to \code{pall$X};
  \code{addpall.GP} also adds \code{Z$y} to \code{pall$Y},
  \code{addpall.CGP} also adds \code{Z$c} to \code{pall$Y},
  and \code{addpall.ConstGP} does both
}
\value{
  nothing is returned, but global variables are modified
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
  \code{\link{PL}}
}
\examples{
## See the demos via demo(package="plgp") and the examples
## section of ?plgp
}
\keyword{ utilities }
