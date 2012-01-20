\name{PL}
\alias{PL}
\alias{plgp}
\title{
Particle Learning Skeleton Method
}
\description{
  Implements the Particle Learning sequential Monte Carlo
  algorithm on the data sequence provided, using re-sample and
  propagate steps
}
\usage{
PL(data, start, end, init, lpredprob, propagate, prior = NULL,
   addpall = NULL, params = NULL, save = NULL, P = 100,
   progress = 10, cont = FALSE, verb = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    function generating the data; for examples see \code{\link{data.GP}}
  }
  \item{start}{
    a scalar \code{integer} specifying the starting \dQuote{time};
    the data entry/sample where PL will start
  }
  \item{end}{
    a scalar \code{integer} specifying the ending \dQuote{time};
    the data entry/sample where PL will stop
  }
  \item{init}{
    function used to initialize the particles at the start of PL;
    for examples see \code{\link{draw.GP}}
  }
  \item{lpredprob}{
    function used to calculate the predictive probability of an
    observation (usually the next one in \dQuote{time}) given a
    particle.  This is the primary function used in the PL re-sample
    step; for examples see \code{\link{lpredprob.GP}}
  }
  \item{propagate}{
    function used to propagate particles given an observation (usually
    the next one in \dQuote{time}); for examples see
    \code{\link{propagate.GP}}
  }
  \item{prior}{
    function used to generate prior parameters that may be
    passed into the \code{data}, \code{init},
    \code{lpredprob} and \code{propagate}
    functions as needed; for examples see \code{\link{prior.GP}}
  }
  \item{addpall}{
    an optional function that adds the new observation (usually
    the next one in \dQuote{time}) to the global \code{pall} variable
    which stores the sufficient information shared by all particles;
    for examples see \code{\link{addpall.GP}}
  }
  \item{params}{
    an optional function called each \code{progress} rounds
    that collects parameters from the particles for
    summary and visualization; for examples see \code{params.GP}
  }
  \item{save}{
    an option function that is called every round to save some
    information about the particles
  }
  \item{P}{
    number of particles to use
  }
  \item{progress}{
    number of PL rounds after which to collect \code{params} and
    draws histograms; a non-positive value or \code{params = NULL}
    skips the progress meter
  }
  \item{cont}{
    if \code{TRUE} then PL will try to use the existing set of particles
    to \dQuote{continue} where it left off; \code{start} and \code{end}
    should be specified appropriately when continuing
  }
  \item{verb}{
    if nonzero, then screen prints will indicate the proportion of PL
    updates finished so far; \code{verb = 1} will cause PL to pause on
    \code{progress} drawings for inspection
  }
}
\details{
  Uses the PL SMC algorithm via the functions provided.  This function
  is just a skeleton framework.  The hard work is in specifying the
  arguments/functions which execute the calculations needed in the
  re-sample and propagate steps.

  PL and uses the global variables \code{pall}, containing
  sufficient information common to all particles, \code{peach},
  containing sufficient information particular to each of the \code{P}
  particles, and \code{psave} containing any saved information.

  Note that PL is designed to be fast for sequential updating
  (of GPs) when new data arrive.  This facilitates efficient sequential
  design of experiments by active learning techniques, e.g.,
  optimization by expected improvement and sequential exploration of
  classification label boundaries by the predictive entropy.  PL is not
  optimized for static inference when all of the data arrive at once,
  in batch
}
\value{
  PL only returns the \code{peach} global variable, containing sufficient
  information particular to each (of the \code{P}) particles
}
\references{
  Carvalho, C., Johannes, M., Lopes, H., and Polson, N. (2008).
  \dQuote{Particle Learning and Smoothing.}
  Discussion Paper 2008-32, Duke University Dept. of Statistical
  Science.
  
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
  
  \url{http://faculty.chicagobooth.edu/robert.gramacy/plgp.html}
}
\author{
  Robert B. Gramacy, \email{rbgramacy@chicagobooth.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{papply}}, \code{\link{draw.GP}},
  \code{\link{data.GP}}, \code{\link{lpredprob.GP}},
  \code{\link{propagate.GP}}, \code{\link{params.GP}},
  \code{\link{pred.GP}}
}
\examples{
## See the demos via demo(package="plgp"); it is important to
## run them with the ask=FALSE argument so that the
## automatically generated plots may refresh automatically
## (without requiring the user to press RETURN)
\dontrun{
## Illustrates regression GPs on a simple 1-d sinusoidal
## data generating mechanism
demo("plgp_sin1d", ask=FALSE)

## Illustrates classification GPs on a simple 2-d exponential
## data generating mechanism
demo("plcgp_exp", ask=FALSE)

## Illustrates classification GPs on Ripley's Cushings data
demo("plcgp_cush", ask=FALSE)

## Illustrates active learning via the expected improvement
## statistic on a simple 1-d data generating mechanism
demo("plgp_exp_ei", ask=FALSE)

## Illustrates active learning via entropy with classification
## GPs on a simple 2-d exponential data generating mechanism
demo("plcgp_exp_entropy", ask=FALSE)

## Illustrates active learning via the integrated expected
## conditional improvement statistic for optimization
## under known constraints on a simple 1-d data generating
## mechanism
demo("plgp_1d_ieci", ask=FALSE)

## Illustrates active learning via the integrated expected
## conditional improvement statistic for optimization under
## unknown constraints on a simple 1-d data generating
## mechanism
demo("plconstgp_1d_ieci", ask=FALSE)

## Illustrates active learning via the integrated expected
## conditional improvement statistic for optimization under
## unknokn constraints on a simple 2-d data generating
## mechanism
demo("plconstgp_2d_ieci", ask=FALSE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{iterations}
