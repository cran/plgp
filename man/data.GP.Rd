\name{data.GP}
\alias{data.GP}
\alias{data.GP.improv}
\alias{data.CGP}
\alias{data.CGP.adapt}
\alias{data.ConstGP}
\alias{data.ConstGP.improv}
\title{
  Supply GP data to PL
}
\description{
  Functions to supply data to PL for Gaussian process (GP)
  regression, classification, or combined unknown constraint
  models
}
\usage{
data.GP(begin, end = NULL, X, Y)
data.GP.improv(begin, end = NULL, f, rect, prior,
               adapt = ei.adapt, cands = 40,
               save = TRUE, oracle = TRUE, verb = 2,
               interp = interp.loess)
data.CGP(begin, end = NULL, X, C)
data.CGP.adapt(begin, end = NULL, f, rect, prior,
               cands = 40, verb = 2, interp=interp.loess)
data.ConstGP(begin, end = NULL, X, Y, C)
data.ConstGP.improv(begin, end = NULL, f, rect, prior,
                    adapt = ieci.const.adapt , cands = 40, 
                    save = TRUE, oracle = TRUE, verb = 2,
                    interp = interp.loess)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{begin}{
    positive \code{integer} starting time for data to be returned
  }
  \item{end}{
    positive \code{integer} (\code{end >= begin}) ending time
    for data being returned; may be \code{NULL} if only data
    at time \code{begin} is needed
  }
  \item{X}{
    \code{data.frame} with at least \code{end} rows containing
    covariates
  }
  \item{Y}{
    vector of length at least \code{end} containing real-valued
    responses
  }
  \item{C}{
    vector of length at least \code{end} containing class labels
  }
  \item{f}{
    function returning a responses when called as \code{f(X)}
    for \code{matrix} \code{X}; for \code{data.GP.improv} the responses
    must be real-valued returned as a vector;
    for \code{data.CGP.adapt} they must be class
    labels returned as a vector;
    for \code{data.ConstGP.improv} they must be pairs of real-valued
    and in \{0,1\} (1 indicates constraint violation), returned as
    a 2-column \code{data.frame}
  }
  \item{rect}{
    bounding rectangle for the inputs \code{X} to \code{f(X)} with
    two columns and rows equalling \code{nrow(X)}
  }
  \item{prior}{
    prior parameters passed from \code{\link{PL}} generated by one of
    the prior functions, e.g., \code{\link{prior.GP}}
  }
  \item{adapt}{
    function that evaluates a sequential design criterion on 
    some candidate locations; the default \code{ei.adapt} EI
    about the minimum; \code{ieci.adapt} providing IECI is another
    possibility , which is hard coded into \code{data.ConstGP.adapt}
  }
  \item{cands}{
    number of Latin Hypercube candidate locations used to choose the
    next adaptively sampled input design point
  }
  \item{save}{
    scalar \code{logical} indicating if the improvment information for
    chosen candidate should be saved in the \code{psave} global variable
  }
  \item{oracle}{
  scalar \code{logical} indicating if the candidates should be
  augmented with the point found to maximize the predictive surface
  (with a search starting at the most recently chosen input)
 }
  \item{verb}{
    verbosity level for printing the progress of improv and other
    adaptive sampling calculations
  }
 \item{interp}{
    function for smoothing of 2-d image plots.  The default comes
    from \code{\link[tgp]{interp.loess}}, but what works best is
    \code{\link[interp]{interp}} which requires the \pkg{interp} or
    \pkg{akima} package
 }
}
\details{
  These functions provide data to PL for Gaussian progress regression
  and classification methods in a variety of ways.  The simplest,
  \code{data.GP} and \code{data.CGP} supply pre-recorded regression and
  classification data stored in data frames and vectors;
  \code{data.ConstGP} is a hybrid that does joint regression and
  classification. The other
  functions provide data by active learning/sequential design:

  The \code{data.GP.improv} function uses expected improvement (EI);
  \code{data.CGP.improv} uses predictive entropy;
  \code{data.ConstGP.improv} 
  uses integrated expected conditional improvement (IECI).  In these
  cases, once the \code{x}-location(s) is/are chosen,
  the function \code{f} is used to provide the response(s)
}
\value{
  The output are vectors or \code{data.frame}s.
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{datagen}
\keyword{methods}
