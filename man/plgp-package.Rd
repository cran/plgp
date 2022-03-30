\name{plgp-package}
\alias{plgp-package}
\docType{package}
\title{Particle Learning of Gaussian Processes}
\description{Sequential Monte Carlo inference for fully Bayesian
  Gaussian process (GP) regression and classification models by
  particle learning (PL).  The sequential nature of inference
  and the active learning (AL) hooks provided facilitate thrifty 
  sequential design (by entropy) and optimization
  (by improvement) for classification and
  regression models, respectively.
  This package essentially provides a generic
  PL interface, and functions (arguments to the interface) which
  implement the GP models and AL heuristics.  Functions for 
  a special, linked, regression/classification GP model and 
  an integrated expected conditional improvement (IECI) statistic 
  is provides for optimization in the presence of unknown constraints.
  Separable and isotropic Gaussian, and single-index correlation
  functions are supported.
  See the examples section of ?plgp and demo(package="plgp") 
  for an index of demos}
\details{
For a fuller overview including a complete list of functions, and
demos, please use \code{help(package="plgp")}.
}
\author{
Robert B. Gramacy \email{rbg@vt.edu}
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
  
  Carvalho, C., Johannes, M., Lopes, H., and Polson, N. (2008).
  \dQuote{Particle Learning and Smoothing}.
  Discussion Paper 2008-32, Duke University Dept. of Statistical
  Science.

  Gramacy, R. (2020).
  \dQuote{Surrogates: Gaussian Process Modeling, Design and Optimization for the Applied Sciences}.
  Chapman Hall/CRC; \url{https://bobby.gramacy.com/surrogates/}

  \url{https://bobby.gramacy.com/r_packages/plgp/}
}
\keyword{ package }
\seealso{ \code{\link{PL}}, \pkg{tgp} }
