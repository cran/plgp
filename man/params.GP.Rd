\name{params.GP}
\alias{params.GP}
\alias{params.CGP}
\alias{params.ConstGP}
\title{
  Extract parameters from GP particles
}
\description{
  Extract parameters from particles for Gaussian process (GP)
  regression, classification, or combined unknown constraint
  models
}
\usage{
params.GP()
params.CGP()
params.ConstGP()
}
%- maybe also 'usage' for other objects documented here.
\details{
  Collects the parameters from each of the particles (contained in
  the global variable \code{peach}) into a
  \code{\link{data.frame}} that can be used for quick
  \code{\link{summary}} and visualization, e.g., via
  \code{\link{hist}}.  These functions are also called to make
  \code{progress} visualizations in \code{\link{PL}}
}
\value{
  returns a \code{data.frame} containing summaries for each
  parameter in its columns
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
Robert B. Gramacy, \email{https://bobby.gramacy.com/r_packages/plgp/}
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
