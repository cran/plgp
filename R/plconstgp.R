#******************************************************************************* 
#
# Particle Learning of Gaussian Processes
# Copyright (C) 2010, University of Cambridge
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
#
#*******************************************************************************


## lpredprob.ConstGP:
##
## for the PL resample step -- combining lpredprob for
## CGP and GP components

lpredprob.ConstGP <- function(z, Zt, prior)
  {
    ## error for hidden constraints
    if(is.na(z$y)) stop("hidden constraints not handled yet")
    
    ## first calculate the GP part as long as z$y is real-valued
    lp <- lpredprob.GP(z, Zt[[1]], prior$GP)

    ## then calculate the CGP part
    lp <- lp + log(pred.CGP(z$x, Zt, prior$CGP, mcreps=1000, cs=z$c))

    ## checks and return
    if(!is.finite(lp)) stop("bad weight")    
    return(lp)
  }


## propagate.ConstGP:
##
## for the PL propagate step -- combining propagate for
## CGP and GP components

propagate.ConstGP <- function(z, Zt, prior)
  {
    ## error for hidden constraints
    if(is.na(z$y)) stop("hidden constraints not handled yet")
    
    ## propagate GP normally if z$y is real-valued, otherwise draw
    Zt[[1]] <- propagate.GP(z, Zt[[1]], prior$GP)

    ## use draw if hidden constraints -- not implemented yet
    ### Zt[[1]] <- draw.GP(Zt[[1]], l=3, h=4, thin=1)

    ## propagate CGP
    Zt <- propagate.CGP(z, Zt, prior$CGP)

    ## return the propagated particle
    return(Zt)
  }


## prior.ConstGP:
##
## default prior specifications for the CGP model and
## the GP model

prior.ConstGP <- function(m, cov.GP=c("isotropic", "separable", "sim"), cov.CGP=cov.GP)
  {
    prior <- list(GP=prior.GP(m, cov.GP), CGP=prior.CGP(m, cov.CGP))
    return(prior)
  }


## draw.ConstGP:
##
## combining draw for CGP and GP components

draw.ConstGP <- function(Zt, prior, l=3, h=4, thin=10)
  {
    ## check if init instead
    if(is.null(Zt)) return(init.ConstGP(prior))

    ## draw for the GP part
    Zt[[1]] <- draw.GP(Zt[[1]], prior$GP, l=l, h=h, thin=thin)

    ## draw for the CGP part
    Zt <- draw.CGP(Zt, prior$CGP, l=l, h=h, thin=thin)

    ## return the fully drawn particle
    return(Zt)
  }


## init.CconstGP:
##
## create a new particle by combining the init functions
## for the CGP and GP

init.ConstGP <- function(prior) 
  {
    Zt <- init.CGP(prior$CGP)
    Zt[[1]] <- init.GP(prior$GP)
    return(Zt)
  }


## pred.ConstGP:
##
## combining the predict functions of GP and CGP

pred.ConstGP <- function(XX, Zt, prior, quants=TRUE)
  {
    GPp <- pred.GP(XX, Zt[[1]], prior$GP, Y=pall$Y, quants=quants)
    CGPp <- pred.CGP(XX, Zt, prior$CGP)
    return(cbind(GPp, CGPp))
  }


## ieci.ConstGP:
##
## combining the predict functions of GP and CGP

ieci.ConstGP <- function(Xcand, Zt, prior, Y=NULL, verb=1)
  {
    ## calculate the predictive probably of the constraint
    ## violation at Xcand under the CGP
    pc2 <- pred.CGP(XX=Xcand, Zt=Zt, prior=prior$CGP, cs=2)

    ## and then caluclate ieci adjusting for pc2
    ieci <- ieci.GP(Xcand=Xcand, Xref=Xcand, Zt[[1]], prior$GP, Y, w=pc2, verb)

    ## calculate in the EI part
    outp <- pred.GP(XX=Xcand, Zt=Zt[[1]], prior=prior$GP, Y=pall$Y, quants=FALSE)

    ## add in the EI part
    ieci <- mean(calc.eis(outp, min(outp$m), pc2)) - ieci

    ## zero out any negative ones (could be -Inf if numerical problems in C)
    ieci[ieci < 0] <- 0

    ## sanity check
    ## if(any(is.nan(ieci))) stop("NaN in ieci")
    
    ## return the EI adjusted IECI
    return(ieci)
  }


## alc.ConstGP:
##
## combining the predict functions of GP and CGP

alc.ConstGP <- function(Xcand, Zt, prior, Y=NULL, verb=1)
  {
    ## calculate the predictive probably of the constraint
    ## violation at Xcand under the CGP
    pc2 <- pred.CGP(XX=Xcand, Zt=Zt, prior=prior$CGP, cs=2)

    ## and then caluclate ieci adjusting for pc2
    alc <- alc.GP(Xcand=Xcand, Xref=Xcand, Zt[[1]], prior$GP, Y, w=pc2, verb)

    ## calculate in the EI part
    outp <- pred.GP(XX=Xcand, Zt=Zt[[1]], prior=prior$GP, Y=pall$Y, quants=FALSE)

    ## add in the ALC part
    alc <- mean(calc.vars(outp, pc2)) - alc

    ## zero out any negative ones (could be -Inf if numerical problems in C)
    alc[alc < 0] <- 0

    ## return the VAR adjusted ALC
    return(alc)
  }


## params.ConstGP:
##
## extracts the params from each GP and each CGP

params.ConstGP <- function()
  {
    ## extract dimensions
    P <- length(peach)
    numGP <- length(peach[[1]])

    ## allocate data frame (DF) to hold parameters
    params <- data.frame(matrix(NA, nrow=P, ncol=3*numGP))

    ## get the names of the parameters, and set them in the DF
    nam <- c()
    for(i in 1:length(peach[[1]])) {
      nam <- c(nam, paste(c("d.", "g.", "lpost."), i, sep=""))
    }
    names(params) <- nam

    ## collect the parameters from the particles
    for(p in 1:P) {
      for(i in 1:numGP) {
        params[p,(i-1)*3+1] <- mean(peach[[p]][[i]]$d)
        params[p,(i-1)*3+2] <- peach[[p]][[i]]$g
        params[p,(i-1)*3+3] <- peach[[p]][[i]]$lpost
      }
    }

    ## return the particles
    return(params)
  }


## data.ConstGP:
##
## extract the appropriate columns from the X matrix, Y vector
## and C vector -- designed to be generic for other cases
## where we would want to get the next observation (end=NULL)
## or a range of observations from begin to end

data.ConstGP <- function(begin, end=NULL, X, Y, C)
  {
    if(is.null(end) || begin == end)
      return(list(x=X[begin,], c=C[begin], y=Y[begin]))
    else if(begin > end) stop("must have begin <= end")
    else return(list(x=as.matrix(X[begin:end,]), y=Y[begin:end], c=C[begin:end]))
  }


## addpall.CGP:
##
## add data to the pall data structure used as utility
## by all particles, combining GP and CGP routines

addpall.ConstGP <- function(Z)
  {
    pall$X <<- rbind(pall$X, Z$x)
    pall$C <<- c(pall$C, Z$c)
    pall$Y <<- c(pall$Y, Z$y)
    pall$D <<- NULL
  }


## data.ConstGP.improv:
##
## use the current state of the particales to calculate
## the next adaptive sample from the posterior predictive
## distribution based on the integrated expected conditional
## improvement (IECI) and the probability that the point
## satisfies the constraint

data.ConstGP.improv <- function(begin, end=NULL, f, rect, prior,
                                adapt=ieci.const.adapt, cands=40, save=TRUE, 
                                oracle=TRUE, verb=2)
  {
    if(!is.null(end) && begin > end) stop("must have begin <= end")
    else if(is.null(end) || begin == end) { ## adaptive sample

      ## choose some adaptive sampling candidates
      Xcand <- lhs(cands, rect)

      ## add a cleverly chosen candidate
      if(oracle) {
        xstars <- findmin.ConstGP(pall$X[nrow(pall$X),], prior)
        xstar <- drop(rectunscale(rbind(xstars), rect))
        Xcand <- rbind(Xcand, xstar)
      }
      
      ## calculate the index with the best IECI
      as <- adapt(Xcand, rect, prior, verb)
      indx <- which.max(as)
      
      ## return the new adaptive sample
      x <- matrix(Xcand[indx,], nrow=1)
      xs <- rectscale(x, rect)

      ## maybe plot something
      if(verb > 1) {
        par(mfrow=c(1,1))
        if(ncol(Xcand) > 1) { ## 2-d+ data
          image(interp(Xcand[,1], Xcand[,2], as))
          points(rectunscale(pall$X, rect))
          points(Xcand, pch=18)
          if(oracle) points(xstar[1], xstar[2], pch=17, col="blue")
          points(x[,1], x[,2], pch=18, col="green")
        } else { ## 1-d data
          o <- order(drop(Xcand))
          plot(drop(Xcand[o,]), as[o], type="l", lwd=2,
               xlab="x", ylab="constrained IECI")
          points(drop(rectunscale(pall$X, rect)), rep(min(as), nrow(pall$X)))
          points(x, min(as), pch=18, col="green")
          if(oracle) points(xstar, min(as), pch=17, col="blue")
          legend("topright", c("chosen point", "oracle candidate"),
                 pch=c(18,17), col=c("green", "blue"), bty="n")
        }
      }

      ## maybe save the max log IECI and xstar
      if(save) {
        if(oracle) psave$xstar <<- rbind(psave$xstar, xstar)
        psave$max.as <<- c(psave$max.as, max(as))
      }

      ## return the adaptively chosen location
      yc <- f(x)
      return(list(x=xs, y=yc$y, c=yc$c))

    } else {  ## create an initial design

      ## calculate a LHS 
      if(verb > 0) cat("initializing with size", end-begin+1, "LHS\n")
      X <- lhs(end-begin+1, rect)
      
      ## get the class labels
      YC <- f(X)
      return(list(x=rectscale(X, rect), y=YC$y, c=YC$c))
    }
  }


## icei.const.adapt:
##
## return the index into Xcand that has the most potential to
## improve the estimate of the minimum via integrated expected
## conitional with constraint probabilities

ieci.const.adapt <- function(Xcand, rect, prior, verb)
  {
    ## calculate the average maximum entropy point
    if(verb > 0)
      cat("taking design point ", nrow(pall$X)+1, " by constained IECI\n", sep="")

    ## adjust the candidates (X) and reference locations (Xref)
    Xcands <- rectscale(Xcand, rect)

    ## get predictive distribution information
    iecis <- papply(Xcand=Xcands, fun=ieci.ConstGP, prior=prior,
                             verb=verb, pre="   IECI")

    ## gather the entropy info for each x averaged over
    ieci <- rep(0, nrow(Xcands))
    for(p in 1:length(iecis)) ieci <- ieci + iecis[[p]]

    ## return the candidate with the most potential
    ## ieci is negated inside ieci.ConstGP
    return(ieci/length(iecis))
  }


## alc.const.adapt:
##
## return the index into Xcand that has the most potential to
## improve reduce the variance via ALC with constraint probabilities

alc.const.adapt <- function(Xcand, rect, prior, verb)
  {
    ## calculate the average maximum entropy point
    if(verb > 0)
      cat("taking design point ", nrow(pall$X)+1, " by constained ALC\n", sep="")

    ## adjust the candidates (X) and reference locations (Xref)
    Xcands <- rectscale(Xcand, rect)

    ## get predictive distribution information
    alcs <- papply(Xcand=Xcands, fun=alc.ConstGP, prior=prior,
                             verb=verb, pre="   ALC")

    ## gather the entropy info for each x averaged over
    alc <- rep(0, nrow(Xcands))
    for(p in 1:length(alcs)) alc <- alc + alcs[[p]]

    ## return the candidate with the most potential
    ## ALC is negated inside alc.ConstGP
    return(alc/length(alcs))
  }


## findmin.ConstGP:
##
## find the minimum of the predictive surface for the MAP particle

findmin.ConstGP <- function(xstart, prior)
  {
    m <- ncol(pall$X)

    ## calculate the MAP particule
    mi <- 1
    for(p in 2:length(peach))
      if(peach[[p]][[1]]$lpost > peach[[mi]][[1]]$lpost) mi <- p
    Zt <- peach[[mi]][[1]]
    
    ## utility for calculations below
    util <- util.GP(Zt, prior$GP, pall$Y, retKi=TRUE)

    ## call the optim function
    xstar <- optim(xstart, pred.mean.GP, method="L-BFGS-B",
                   lower=rep(0,m), upper=rep(1,m), util=util,
                   cov=util$GP$cov, dparam=Zt$d, gparam=Zt$g)$par
    return(xstar)
  }
      
