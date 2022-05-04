#######################################################################
## Copyright 2008, Novartis Pharma AG
##
## This program is Open Source Software: you can redistribute it
## and/or modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see http://www.gnu.org/licenses/.

# functions to design trials for MCPMod
library(mvtnorm)
library(lattice)
library(shiny)
guesst <-
  function(d, p, model = c("emax", "exponential", "logistic", "quadratic",
                           "betaMod", "sigEmax"), less = TRUE,  local = FALSE, dMax, Maxd, scal)
  {
    ## get guesstimates for standardized model parameters(s)
    ## d = dose corresponding to prior expected percentage "p" of maximum effect
    ## p = percentage corresponding to "d"
    ## model = name of dose response model
    ## less = logical only needed in case of quadratic model.
    ##        determines if d is smaller (less=T) or larger (less=F)
    ##        than d_opt
    ## local = logical indicating whether local (ie exact) or asymptotic versions
    ##         of guesstimate should be derived for emax, logistic and
    ##         sigEmax model (for other models method is exact)
    ##         Branson et al. derive the asymptotic
    ##         expressions (ie. assuming that f^0(Maxd)=1)
    ##         in case local version is used, convergence problems
    ##         with underlying nonlinear optimization can occur
    ## dMax = dose at which max effect occurs (needed only for betaMod)
    ## Maxd = Maximum dose administered
    ## scal = dose scale (needed only for betaMod)
    model <- match.arg(model)
    if(any(p <= 0) | any(p > 1)) stop("must have 0 < p <= 1")
    switch(model,
           emax = {
             if(!local){
               c(ed50 = mean(d * (1 - p)/p))
             } else {
               if (any(p <= d/Maxd))
                 stop("must have p > d/Maxd, for local version")
               val <- (d/p-d)/(1-d/(Maxd*p))
               c(ed50=mean(val))
             }
           },
           exponential = {
             if(any(p >= d/Maxd)) stop("must have p < d/Maxd")
             init <- d/log(1 + p)
             fooexp <- function(delta, d, p, Maxd){
               sum((exponential(d, 0, 1, delta)/
                      exponential(Maxd, 0, 1, delta) - p)^2)
             }
             val <- optimize(fooexp, c(0, 2*Maxd), d=d, p=p, Maxd=Maxd)$minimum
             c(delta = mean(val))
           },
           logistic = {
             if(length(d) == 1) {
               stop("logistic model needs at least two pairs (d,p)")
             }
             logit <- function(p) log(p/(1-p))
             if(length(d) == 2) {
               ed50  <- diff(rev(d)*logit(p))/diff(logit(p))
               delta <- diff(d)/diff(logit(p))
               res <- c(ed50 = ed50, delta = delta)
             } else {
               m <- lm(logit(p)~d)
               par <- coef(m)
               names(par) <- NULL
               res <- c(ed50 = -par[1]/par[2], delta = 1/par[2])
             }
             if(local){
               foolog <- function(par, d, p, Maxd) {
                 e0 <- logistic(0,0,1,par[1],par[2])
                 sum(((logistic(d,0,1,par[1],par[2]) - e0)/
                        (logistic(Maxd,0,1,par[1],par[2])-e0)-p)^2)
               }
               res <- try(optim(par=res, fn=foolog, d=d, p=p, Maxd=Maxd))
               if(res$convergence > 0)
                 stop("cannot find guesstimates for specified values")
               else res <- res$par
             }
             if(res[1] < 0)
               message("specified values lead to negative ed50, which should be positive")
             res
           },
           quadratic = {
             aux <- sqrt(1 - p)
             if (less) c(delta = mean(-(1 - aux)/(2 * d)))
             else  c(delta = mean(-(1 + aux)/(2 * d)))
           },
           betaMod = {
             if(scal <= dMax)
               stop("scal needs to be larger than dMax to calculate guesstimate")
             k <- dMax/(scal-dMax)
             val <- d^k*(scal-d)/(dMax^k*(scal-dMax))
             beta <- log(p)/(log(val))
             c(delta1 = mean(k*beta), delta2 = mean(beta))
           },
           sigEmax = {
             if(length(d) == 1) {
               stop("sigmoid Emax model needs at least two pairs (d,p)")
             }
             if(length(d) == 2){
               num <- log((p[1]*(1-p[2]))/(p[2]*(1-p[1])))
               h <- num/log(d[1]/d[2])
               ed50 <- ((1-p[1])/p[1])^(1/h)*d[1]
               res <- c(ed50=ed50, h=h)
             } else {
               y <- log((1-p)/p)
               x <- log(d)
               par <- coef(lm(y~x))
               names(par) <- NULL
               res <- c(ed50 = exp(par[1]/-par[2]), delta = -par[2])
             }
             if(local) {
               fooSE <- function(par, d, p, Maxd) {
                 sum((sigEmax(d,0,1,par[1],par[2])/
                        sigEmax(Maxd,0,1,par[1],par[2])-p)^2)
               }
               res <- try(optim(par=res, fn=fooSE, d=d, p=p, Maxd=Maxd))
               if(res$convergence > 0)
                 stop("cannot find guesstimates for specified values")
               else res <- res$par
             }
             if(res[1] < 0)
               message("specified values lead to negative ed50, which should be positive")
             res
           })
  }
### example call
### guesst(0.5, 0.8, "emax")
### guesst(0.5, 0.8, "emax", Maxd=1, local=T)
### guesst(0.5, 0.4, "exponential", Maxd = 1)
### guesst(0.7, 1, "quadratic") # maximum response at dose 0.7
### guesst(c(0.4, 0.7), c(0.5, 0.95), "logistic")
### guesst(c(0.4, 0.7), c(0.5, 0.95), "logistic", local = T, Maxd = 1)
### guesst(c(0.4, 0.7), c(0.5, 0.95), "sigEmax")
### guesst(c(0.4, 0.7), c(0.5, 0.95), "sigEmax", local = T, Maxd = 1)
### guesst(0.4,0.5, "betaMod", dMax = 0.8, scal=1.3)

linear <-
  function(dose, e0, delta) e0 + delta * dose

linlog <-
  function(dose, e0, delta, off = 1) linear(log(dose + off), e0, delta)

emax <-  function(dose, e0, eMax, ed50){
  sigEmax(dose, e0, eMax, ed50, 1)
}

quadratic <-
  function(dose, e0, b1, b2) e0 + b1 * dose + b2 * dose^2

exponential <- function(dose, e0, e1, delta){
  e0 + e1*(exp(dose/delta) - 1)
}

logistic <-
  function(dose, e0, eMax, ed50, delta)
  {
    e0 + eMax/(1 + exp((ed50 - dose)/delta))
  }

betaMod <-
  function(dose, e0, eMax, delta1, delta2, scal)
  {
    maxDens <- (delta1^delta1)*(delta2^delta2)/
      ((delta1 + delta2)^(delta1+delta2))
    dose <- dose/scal
    e0 + eMax/maxDens * (dose^delta1) * (1 - dose)^delta2
  }

sigEmax <-
  function(dose, e0, eMax, ed50, h)
  {
    e0 + eMax*dose^h/(ed50^h + dose^h)
  }

## means corresponding to different models and doses
modelMeans <-
  ## generate mean vectors for models and initial values defined by user
  function(models, doses, std = TRUE, off = 0.1*max(doses), scal = 1.2*max(doses))
  {
    ## models - list with names given by the model functions to be used
    ##          and components as eitheir vectors or matrices with the
    ##          parameter values ; assumption is that first argument to
    ##          model function is always the doses
    ## doses  - doses to be used in design
    ## std    - logical variable indicating whether to assume
    ##          standardized versions of built-in models

    namMod <- names(models)
    if(length(namMod) != length(unique(namMod))){
      stop("only one list entry allowed for each model class in 'models' argument")
    }
    ## built-in models with methods for standardized versions
    biModels <- c("emax", "linlog", "linear", "quadratic",
                  "exponential", "logistic", "betaMod", "sigEmax")
    nModels <- length(models)             # number of model elements
    contMap <- vector("list", nModels)
    names(contMap) <- namMod
    val <- list()
    k <- 1
    nams <- character()
    for(nm in namMod) {
      pars <- models[[nm]]
      if (!is.null(pars) && !is.numeric(pars)) {
        stop("elements of \"models\" must be NULL or numeric")
      }
      if (is.element(nm, biModels) && std) {
        if (!is.element(nm, c("betaMod","logistic","sigEmax"))) {
          ## all others have single std parameter
          if (is.matrix(pars)) {
            stop("For standardized ", nm,
                 " models element cannot be matrix")
          }
          nmod <- length(pars)
          if (nmod > 1) {                 # multiple models
            ind <- 1:nmod
            nams <- c(nams, paste(nm, ind, sep = ""))
            contMap[[nm]] <- (k - 1) + ind
            for(j in 1:nmod) {
              if (nm == "linlog") pars1 <- c(0, 1, off)
              else pars1 <- c(0, 1, pars[j])
              val[[k]] <- do.call(nm, c(list(doses), as.list(pars1)))
              k <- k + 1
            }
          } else {                        # single model
            nams <- c(nams, nm)
            contMap[[nm]] <- k
            if (nm == "quadratic") names(pars) <- NULL
            if (nm == "linlog") pars <- c(0, 1, off)
            else pars <- c(0, 1, pars)
            val[[k]] <- do.call(nm, c(list(doses), as.list(pars)))
            k <- k + 1
          }
        } else {                          # logistic, betaMod
          if (is.matrix(pars)) {          # multiple models
            if (ncol(pars) != 2) {
              stop("standardized ", nm," model must have two parameters")
            }
            nmod <- nrow(pars)            # number of models
            ind <- 1:nmod
            nams <- c(nams, paste(nm, ind, sep = ""))
            contMap[[nm]] <- (k - 1) + ind
            for(j in 1:nmod) {
              if(nm == "betaMod")     pars1 <- c(0, 1, pars[j, ], scal)
              else pars1 <- c(0, 1, pars[j, ])
              val[[k]] <- do.call(nm, c(list(doses), as.list(pars1)))
              k <- k + 1
            }
          } else {                        # single model
            if (length(pars) != 2) {
              stop("standardized ", nm," model must have two parameters")
            }
            nams <- c(nams, nm)
            contMap[[nm]] <- k
            if(nm == "betaMod") pars <- c(pars, scal)
            val[[k]] <-
              do.call(nm, c(list(doses), as.list(c(0, 1, pars))))
            k <- k + 1
          }
        }
      } else {                        # user-defined or non-standardized
        if (is.matrix(pars)) {            # multiple models
          nmod <- nrow(pars)              # number of models
          if(nm == "linlog")  pars <- cbind(pars, off)
          if(nm == "betaMod") pars <- cbind(pars, c(scal))
          ind <- 1:nmod
          nams <- c(nams, paste(nm, ind, sep = ""))
          contMap[[nm]] <- (k - 1) + ind
          for(j in 1:nmod) {
            val[[k]] <- do.call(nm, c(list(doses), as.list(pars[j,])))
            k <- k + 1
          }
        } else {                      # single model
          if(nm == "linlog")  pars <- c(pars, off)
          if(nm == "betaMod") pars <- c(pars, scal)
          nams <- c(nams, nm)
          contMap[[nm]] <- k
          val[[k]] <- do.call(nm, c(list(doses), as.list(pars)))
          k <- k + 1
        }
      }
    }
    muMat <- do.call("cbind", val)
    dimnames(muMat) <- list(doses, nams)
    attr(muMat, "contMap") <- contMap
    muMat
  }

## Model contrasts
getStdCont <-
  function(cont)
  {
    ## center and normalize vector
    cont <- cont - mean(cont)
    cont/sqrt(sum(cont^2))
  }

### simplified version
getOptCont <-
  function(mu, n = 1)
  {
    mn <- sum(mu * n)/sum(n)
    getStdCont(n * (mu - mn))
  }

modContr <-
  function(means, n = 1)
  {
    ## obtains model contrasts corresponding to model means and
    ## sample sizes n
    apply(means, 2, getOptCont, n = n)
  }

# control parameters for functions in package mvtnorm
mvtnorm.control <-
  function(maxpts = 30000, abseps = 0.001,
           interval = NULL, releps = 0)
  {
    out <- list(interval = interval, maxpts = maxpts, abseps = abseps,
                releps = releps)
    class(out) <- "GenzBretz"
    out
  }

critVal <-
  function(cMat, n, alpha = 0.025, control = mvtnorm.control(),
           twoSide = FALSE, corMat = NULL)
  {
    ## function to calculate critical value
    ##
    ## cMat - model contrast matrix nDose x nMod
    ## n - scalar or vector with sample size(s)
    ## alpha - significance level
    ## control - see mvtnorm.control function
    ## twoSide - logical variable indicating if one or two-sided test
    ## corMat - correlation matrix

    nD <- nrow(cMat)
    nMod <- ncol(cMat)
    if (length(n) == 1) {
      nDF <- nD * (n - 1)
    } else {
      if (length(n) != nD) {
        stop("sample size vector must have length equal to number of doses")
      }
      nDF <- sum(n) - nD
    }
    if(is.null(corMat)){
      corMat <- t(cMat)%*%(cMat/n)
      den  <- sqrt(crossprod(t(colSums(cMat^2/n))))
      corMat <- corMat / den
    }
    if(twoSide) tail <- "both.tails"
    else tail <- "lower.tail"
    if (!missing(control)) {
      if(!is.list(control)) {
        stop("when specified, 'control' must be a list")
      }
      ctrl <- do.call("mvtnorm.control", control)
    } else {
      ctrl <- control
    }
    qmvtCall <- c(list(1-alpha, tail = tail, df = nDF, corr = corMat,
                       algorithm=ctrl, interval=ctrl$interval))
    do.call("qmvt", qmvtCall)$quantile
  }

planMM <-
  function(models, doses, n, off = 0.1*max(doses), scal = 1.2*max(doses), std = TRUE,
           alpha = 0.025, twoSide = FALSE, control = mvtnorm.control(),
           cV = TRUE, muMat = NULL)
  {
    ## returns the critical value, the contrast and correlation matrices
    ## models - list with names given by the model functions to be used
    ##          and components as either vectors or matrices with the
    ##          parameter values ; assumption is that first argument to
    ##          model function is always the doses
    ## doses  - doses to be used in design
    ## n - scalar or vector with sample size(s)
    ## off, scal - additional parameters for the linear in log and beta model
    ## std    - logical variable indicating whether to assume
    ##          standardized versions of built-in models
    ## alpha - significance level
    ## twoSide - logical variable indicating if one or two-sided test
    ## control - control arguments for mvtnorm function(s)
    ## cV - logical indicating whether critical value should be calc.
    ## muMat - an optional matrix with means as columns, dimnames should
    ##         be given (dose levels and names of contrasts)

    if (missing(models)) {
      ## may pass necessary information via muMat
      if (is.null(muMat)) {
        stop("either models or muMat need to be specified")
      }
      mu <- muMat
    } else {
      mu <- modelMeans(models, doses, std, off, scal)
    }
    contMat <- modContr(mu, n)
    corMat <- t(contMat)%*%(contMat/n)
    den  <- sqrt(crossprod(t(colSums(contMat^2/n))))
    corMat <- corMat / den
    if(cV){
      critV <- critVal(contMat, n, alpha, control = control, twoSide = twoSide)
    }
    else critV <- NULL
    res <- list(contMat = contMat, critVal = critV, muMat = mu,
                corMat = (corMat+t(corMat))/2 )
    attr(res, "n") <- n
    attr(res, "alpha") <- alpha
    if (twoSide) {
      attr(res, "side") <- "two-sided"
    } else {
      attr(res, "side") <- "one-sided"
    }
    oldClass(res) <- "planMM"
    res
  }

print.planMM <-
  function(x, digits = 3, ...)
  {
    cat("MCPMod planMM\n")
    cat("\n","Optimal Contrasts:","\n", sep="")
    print(round(x$contMat, digits))
    if(!is.null(x$critVal)){
      names(x$critVal) <- NULL
      cat("\n","Critical Value ", sep="")
      cat(paste("(alpha = ", attr(x, "alpha"),", ", attr(x, "side"), "): ",
                round(x$critVal, digits), sep=""),"\n")
    }
    cat("\n","Contrast Correlation Matrix:","\n", sep="")
    print(round(x$corMat, digits))
    cat("\n")
  }

### example call (example from JBS paper)
### doses <- c(0, 10, 25, 50, 100, 150)
### models <- list(linear = NULL, emax = c(25),
###                logistic = c(50, 10.88111), exponential = c(85),
###                betaMod = matrix(c(0.33, 2.31, 1.39, 1.39), byrow=T,nrow=2))
### planMM(models, doses, 50, scal = 200, alpha = 0.05) # compare with JBS 16, 645-646
### example with muMat
### dvec <- c(0, 10, 50, 100)
### mu1 <- c(1, 2, 2, 2)
### mu2 <- c(1, 1, 2, 2)
### mu3 <- c(1, 1, 1, 2)
### mMat <- cbind(mu1, mu2, mu3)
### dimnames(mMat)[[1]] <- dvec
### planMM(muMat = mMat, doses = dvec, n = 30)

plot.planMM <- function (x, superpose = TRUE, xlab = "Dose",
                         ylab = NULL, resp = c("contrasts", "means"), ...)
{
  resp <- match.arg(resp)
  if (is.null(ylab)) {
    if (resp == "contrasts") {
      ylab <- "Contrast coefficients"
    }
    else {
      ylab <- "Normalized model means"
    }
  }
  cM <- x$contMat
  if (resp == "means") {
    cM <- t(t(x$muMat)/apply(x$muMat, 2, max))
  }
  nD <- nrow(cM)
  nM <- ncol(cM)
  cMtr <- data.frame(resp = as.vector(cM), dose = rep(as.numeric(dimnames(cM)[[1]]),
                                                      nM), model = factor(rep(dimnames(cM)[[2]], each = nD),
                                                                          levels = dimnames(cM)[[2]]))
  if (superpose) {
    spL <- trellis.par.get("superpose.line")
    spL$lty <- rep(spL$lty, nM%/%length(spL$lty) + 1)[1:nM]
    spL$lwd <- rep(spL$lwd, nM%/%length(spL$lwd) + 1)[1:nM]
    spL$col <- rep(spL$col, nM%/%length(spL$col) + 1)[1:nM]
    xyplot(resp ~ dose, data = cMtr, subscripts = TRUE, groups = cMtr$model,
           panel = panel.superpose, type = "o", xlab = xlab,
           ylab = ylab, key = list(lines = spL, transparent = TRUE,
                                   text = list(levels(cMtr$model), cex = 0.9), columns = ifelse(nM <
                                                                                                  5, nM, min(4,ceiling(nM/min(ceiling(nM/4),3))))), ...)
  }
  else {
    xyplot(resp ~ dose | model, data = cMtr, type = "o",
           xlab = xlab, ylab = ylab, strip = function(...) strip.default(...,
                                                                         style = 1), ...)
  }
}

### doses <- c(0, 10, 25, 50, 100, 150)
### models <- list(linear = NULL, emax = c(25),
###                logistic = c(50, 10.88111), exponential = c(85),
###                betaMod = matrix(c(0.33, 2.31, 1.39, 1.39), byrow=T,nrow=2))
### pM <- planMM(models, doses, 50, scal = 200)
### plot(pM)
### plot(pM, superpose=F, xlab="Different axis name")
### plot(pM, resp = "means")
### #example with muMat
### dvec <- c(0, 10, 50, 100)
### mu1 <- c(1, 2, 2, 2)
### mu2 <- c(1, 1, 2, 2)
### mu3 <- c(1, 1, 1, 2)
### mMat <- cbind(mu1, mu2, mu3)
### dimnames(mMat)[[1]] <- dvec
### pM <- planMM(muMat = mMat, doses = dvec, n = 30)
### plot(pM)
### plot(pM, superpose=F, xlab="Different axis name")
### plot(pM, resp = "means")

getPars <-
  function(model, doses, initEstim, base, maxEff, off = 0.1*max(doses), scal = 1.2*max(doses))
  {
    ## calculates the location and scale parameters corresponding to
    ## given base, maxEff, and guesstimates
    ## model - a character string with the model name. Built-in models
    ##         have their full parameterization derived internally. For
    ##         user-defined models, it is assumed the existence of a
    ##         function named as "Par" appended to end of model (e.g.,
    ##         for model = "cubic", it is assumed that there is function
    ##         cubicPar that calculates the necessary parameters; this
    ##         function is assumed to have arguments doses, initEstim,
    ##         base, and maxEff, in that order. (see below for an example)
    ## doses - doses to be used in design
    ## base, maxEff - baseline and maximum effect (over placebo)
    ## initEstim - vector of guesstimates
    ## off, scal - additional parameters for the linear in log and beta model
    switch(model,
           linear = {
             e1 <- maxEff/max(doses)
             Par <- c(base, e1)
           },
           linlog = {
             e1 <- maxEff/(log(max(doses) + off) - log(off))
             Par <- c(base-e1*log(off), e1)
           },
           quadratic = {
             dMax <- 1/(-2*initEstim)
             b1 <- maxEff/(dMax + initEstim*dMax^2)
             b2 <- initEstim * b1
             Par <- c(base, b1, b2)
           },
           emax = {
             emax.p <- maxEff * (initEstim + max(doses))/max(doses)
             Par <- c(base, emax.p, initEstim)
           },
           exponential = {
             e1 <- maxEff/(exp(max(doses)/initEstim) - 1)
             e0 <- base
             Par <- c(e0, e1, initEstim)
           },
           logistic = {
             emax.p <- maxEff/
               (logistic(max(doses),0,1, initEstim[1], initEstim[2]) -
                  logistic(0, 0, 1, initEstim[1], initEstim[2]))
             e0 <- base-emax.p*logistic(0,0,1,initEstim[1], initEstim[2])
             Par <- c(e0, emax.p, initEstim[1], initEstim[2])
           },
           betaMod = {
             Par <- c(base, maxEff, initEstim)
           },
           sigEmax = {
             ed50 <- initEstim[1]
             h <- initEstim[2]
             dmax <- max(doses)
             eMax <- maxEff*(ed50^h+dmax^h)/dmax^h
             Par <- c(e0 = base, eMax = eMax, ed50 = ed50, h = h)
           },
           {
             Par <- do.call(paste(model,"Par", sep=""),
                            list(doses, initEstim, base, maxEff))
           }
    )
    Par
  }

### example call
### doses <- c(0, 10, 25, 50, 100, 150)
### getPars("emax", doses, 25, 0, 0.4)
### getPars("logistic", doses, c(50, 10.88111), 0, 0.4) # compare JBS 16, p.650
### getPars("betaMod", doses, initEstim = c(0.33, 2.31), base = 0,
###         maxEff = 0.4)
### #example for user model
### sigEmax2 <- function(dose, e0, eMax, ed50, h){
###   e0 + eMax * ( dose^h / (ed50^h + dose^h) )
### }
### # function to return location and scale parameters
### sigEmax2Par <- function(dose, initEstim, base, maxEff){
###   # function to get linear parameters
###   # ed50 parameter assumed to be first in initEstim
###   ed50 <- initEstim[1]
###   h <- initEstim[2]
###   dmax <- max(dose)
###   emax <- maxEff*(ed50^h+dmax^h)/dmax^h
###   c(base, emax, initEstim)
### }
### getPars("sigEmax2", doses, initEstim = c(50,2), base = 0, maxEff = 1)


fullMod <-
  function(models, doses, base, maxEff, off = 0.1*max(doses), scal = 1.2*max(doses))
  {
    ## calculates parameters for all models in the candidate set
    ## returns a list (similar to the models list) but with all model parameters.
    ## models - list with names given by the model functions to be used
    ##          and components as either vectors or matrices with the
    ##          parameter values
    ## doses  - doses to be used in design
    ## base, maxEff - baseline, maximum effect (over placebo)
    ## off, scal - additional parameters for the linear in log and beta model

    namMod <- names(models)
    if(length(namMod) != length(unique(namMod))){
      stop("only one list entry allowed for each model class in 'models' argument")
    }
    complModels <- list()
    i <- 0
    for(nm in namMod){
      pars <- models[[nm]]
      if(is.null(pars)){
        Pars <- getPars(nm, doses, NULL, base, maxEff, off); i <- i+1
      }
      else{
        if(!is.element(nm,c("logistic", "betaMod", "sigEmax"))){
          nmod <- length(pars)
          if(nmod > 1){
            Pars <- matrix(ncol=3, nrow=nmod)
            for(j in 1:length(pars)){
              Pars[j,] <-  getPars(nm, doses, as.vector(pars[j]), base, maxEff)
            }
            i <- i+1
          }
          else {
            Pars <-  getPars(nm, doses, as.vector(pars), base, maxEff)
            i <- i+1
          }
        }
        else {
          if(is.matrix(pars)){
            Pars <- matrix(ncol=4, nrow=nrow(pars))
            for(j in 1:nrow(pars)){
              Pars[j,] <-  getPars(nm, doses, as.vector(pars[j,]), base, maxEff)
            }
            i <- i+1
          }
          else {
            Pars <-  getPars(nm, doses, as.vector(pars), base, maxEff); i <- i+1
          }
        }
      }
      complModels[[i]] <- Pars
    }
    names(complModels) <- names(models)
    ## store information for use later
    attr(complModels, "doses") <- doses
    attr(complModels, "base") <- base
    attr(complModels, "maxEff") <- maxEff
    attr(complModels, "off") <- off
    attr(complModels, "scal") <- scal
    oldClass(complModels) <- "fullMod"
    complModels
  }

### example call
### doses <- c(0, 10, 25, 50, 100, 150)
### models <- list(linear = NULL, emax = c(25),
###                logistic = c(50, 10.88111), exponential = c(85),
###                betaMod = matrix(c(0.33, 2.31, 1.39, 1.39), byrow=T,nrow=2))
### fullMod(models, doses, base = 0, maxEff = 0.4, scal = 200)

plotModels <- function (models, doses, base, maxEff, nPoints = 200, off = 0.1*max(doses),
                        scal = 1.2 * max(doses), superpose = FALSE, ylab = "Model means",
                        xlab = "Dose", ...)
{
  if (inherits(models, "fullMod")) {
    doses <- attr(models, "doses")
    base <- attr(models, "base")
    maxEff <- attr(models, "maxEff")
    off <- attr(models, "off")
    scal <- attr(models, "scal")
    Models <- models
  }
  else {
    Models <- fullMod(models, doses, base, maxEff, off, scal)
  }
  doseSeq <- sort(union(seq(min(doses), max(doses), length = nPoints),
                        doses))
  resp <- modelMeans(Models, doseSeq, std = FALSE, off, scal)
  nams <- dimnames(resp)[[2]]
  nM <- length(nams)
  if (superpose) {
    respdata <- data.frame(response = c(resp),
                           dose = rep(doseSeq, ncol(resp)),
                           model = factor(rep(nams, each = length(doseSeq)),
                                          levels = nams))
    spL <- trellis.par.get("superpose.line")
    spL$lty <- rep(spL$lty, nM%/%length(spL$lty) + 1)[1:nM]
    spL$lwd <- rep(spL$lwd, nM%/%length(spL$lwd) + 1)[1:nM]
    spL$col <- rep(spL$col, nM%/%length(spL$col) + 1)[1:nM]
    xyplot(response ~ dose, data = respdata, subscripts = TRUE,
           groups = respdata$model, panel.data = list(base = base, maxEff = maxEff,
                                                      doses = doses), xlab = xlab, ylab = ylab, panel = function(x,
                                                                                                                 y, subscripts, groups, ..., panel.data) {
                                                        panel.abline(h = c(panel.data$base, panel.data$base +
                                                                             panel.data$maxEff), lty = 2)
                                                        panel.superpose(x, y, subscripts, groups, type = "l",
                                                                        ...)
                                                        ind <- !is.na(match(x, panel.data$doses))
                                                        panel.superpose(x[ind], y[ind], subscripts[ind],
                                                                        groups, ...)
                                                      }, key = list(lines = spL, transparent = TRUE, text = list(nams,
                                                                                                                 cex = 0.9), columns = ifelse(nM <
                                                                                                                                                5, nM, min(4,ceiling(nM/min(ceiling(nM/4),3))))), ...)
  }
  else {
    respdata <- data.frame(response = c(resp), dose = rep(doseSeq,
                                                          ncol(resp)), model = factor(rep(nams, each = length(doseSeq))))
    xyplot(response ~ dose | model, data = respdata, panel.data = list(base = base,
                                                                       maxEff = maxEff, doses = doses), xlab = xlab, ylab = ylab,
           panel = function(x, y, ..., panel.data) {
             panel.abline(h = c(panel.data$base, panel.data$base +
                                  panel.data$maxEff), lty = 2)
             panel.xyplot(x, y, type = "l", ...)
             ind <- match(panel.data$doses, x)
             panel.xyplot(x[ind], y[ind], ...)
           }, strip = function(...) strip.default(..., style = 1),
           as.table = TRUE,...)
  }
}

plot.fullMod <-
  function(x, ...)
  {
    plotModels(x, ...)
  }

### example call
### doses <- c(0,10,25,50,100,150)
### models <- list(linear = NULL, emax = c(25),
###                logistic = c(50, 10.88111), exponential = c(85),
###                betaMod = matrix(c(0.33, 2.31, 1.39, 1.39), byrow=T,nrow=2))
### plotModels(models, doses, base = 0, maxEff = 0.4, scal = 200)
### plotModels(models, doses, base = 0, maxEff = 0.4, scal = 200, superpose = T)

powCalc <-
  function(cMat, n, alpha = 0.025, delta = NULL, mu = NULL,
           sigma = NULL, cVal = NULL, corMat = NULL,
           twoSide = FALSE, control = mvtnorm.control())
  {
    ## cMat - contrast matrix with rows (as number of doses) and columns
    ##        (as number of models) - exactly as in the output of planMM
    ## n - scalar or vector with sample size(s)
    ## alpha - significance level
    ## delta - non-centrality vector (in same order as models)
    ## mu - vector of model means (either delta or mu and sigma need to be given)
    ## sigma - standard deviation of response (to build delta vector)
    ## cVal - optional critical value for the test
    ## twoSide - logical variable indicating if one or two-sided test

    nD <- nrow(cMat)
    nMod <- ncol(cMat)
    if (length(n) == 1) {
      nDF <- nD * (n - 1)
    } else {
      if (length(n) != nD) {
        stop("sample size vector must have length equal to number of doses")
      }
      nDF <- sum(n) - nD
    }
    if (is.null(delta)) {
      den  <- sqrt(colSums(cMat^2/n))
      delta <- as.vector((t(cMat) %*% mu)/(den*sigma))
    }
    if(is.null(corMat)){
      corMat <- t(cMat)%*%(cMat/n)
      den  <- sqrt(crossprod(t(colSums(cMat^2/n))))
      corMat <- corMat / den
    }
    if (is.null(cVal)) {
      cVal <-
        critVal(cMat, n, alpha, control, twoSide, corMat)
    }
    if(twoSide){
      lower <- rep(-cVal, nMod)
    }
    else lower <- rep(-Inf, nMod)
    upper <- rep(cVal, nMod)
    if (!missing(control)) {
      if(!is.list(control)) {
        stop("when specified, 'control' must be a list")
      }
      ctrl <- do.call("mvtnorm.control", control)
    } else {
      ctrl <- control
    }
    ctrl$interval <- NULL      # not used with pmvt
    pmvtCall <- c(list(lower, upper, df = nDF, corr = corMat,
                       delta = delta, algorithm=ctrl))
    as.vector(1 - do.call("pmvt", pmvtCall))
  }
### example call
### doses <- c(0,10,25,50,100,150)
### models <- list(linear = NULL, emax = c(25),
###                logistic = c(50, 10.88111), exponential=c(85),
###                betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=T, nrow=2))
### plMM <- planMM(models, doses, 50, scal = 200, alpha = 0.05)
### compMod <- fullMod(models, doses, base = 0, maxEff = 0.4, scal = 200)
### muMat <- modelMeans(compMod, doses, F, scal = 200)
### powCalc(plMM$contMat, 50, mu = muMat[,3], sigma = 1, cVal = plMM$critVal, control = list(maxpts=100000)) # Power for logistic model
### powCalc(plMM$contMat, 50, mu = muMat[,1], sigma = 1, cVal = plMM$critVal) # Power for linear model
### # compare with JBS 16, 650

powerMM <-
  function(models, doses, base, maxEff, sigma, lower, upper,
           step, sumFct = c("min", "mean", "max"), off = 0.1*max(doses),
           scal = 1.2*max(doses), alpha = 0.025, twoSide = FALSE,
           control = mvtnorm.control(), muMat = NULL, alRatio = NULL,
           typeN = c("arm", "total"), ...)
  {
    ## generates a matrix with power values for different sample sizes
    ## (balanced case only) and models
    ## models - either a list with model parameters as in planMM
    ##          or an object of class fullMod
    ## doses, maxEff, base, off, scal - as in getPars function
    ## sigma - expected standard deviation
    ## lower, upper - maximum and minimum group sample size for which the
    ##                power is plotted
    ## step - in which steps should the power be calculated
    ##        (it is calculated at seq(lower,upper,by=step))
    ## sumFct - a character vector giving the names of the summary
    ##          functions to be called (each function should return a numeric
    ##          of length 1)
    ## off,scal - additional parameters for linlog/beta model
    ## alpha - level of significance
    ## twoSide - logical variable indicating if one or two-sided test
    ## control - control list for mvtnorm fcts., see mvtnorm.control()
    ## muMat - an optional matrix with means as columns
    ## alRatio - optional vector with sample size allocation ratios
    ##           must all be > 0; lower and upper are interpreted as
    ##           limit values for dose(s) with smallest alloc. ratios
    ## ... - possible additional arguments for sumFct

    if (!missing(models)) {
      ## if models is of class fullMod, will ignore other arguments
      if (inherits(models, "fullMod")) {
        doses <- attr(models, "doses")
        base <- attr(models, "base")
        maxEff <- attr(models, "maxEff")
        off <- attr(models, "off")
        scal <- attr(models, "scal")
        Models <- models
      } else {
        Models <- fullMod(models, doses, base, maxEff, off, scal)
      }
      muMat <- modelMeans(Models, doses, FALSE, off, scal)
    } else {
      if (is.null(muMat)) {
        stop("muMat must be specified, when models is not")
      }
    }
    nD <- length(doses)
    # allocation ratios
    Ntype <- match.arg(typeN)
    if(!is.null(alRatio)){
      if(length(alRatio)!= nrow(muMat)){
        stop("alRatio needs to be of same length as number of dose levels")
      }
      if(any(alRatio <= 0)){
        stop("all entries of alRatio need to be positive")
      }
      if (Ntype == "arm") {
        alRatio <- alRatio/min(alRatio)
      } else {
        alRatio <- alRatio/sum(alRatio)
      }
    } else {
      alRatio <- 1
    }

    contMat <- modContr(muMat, upper*alRatio)
    nmod <- ncol(contMat)
    pow <- numeric(nmod)
    j <- 1
    nSeq <- seq(lower, upper, by = step)
    powMat <- matrix(nrow=length(nSeq), ncol=nmod)
    for(n in nSeq) {
      N <- round(n * alRatio)             # simple rounding
      cVal <- critVal(contMat, N, alpha, control, twoSide)
      for(i in 1:nmod) {
        den  <- sqrt(colSums(contMat^2/N))
        delta <- as.vector((t(contMat) %*% muMat[,i])/(den*sigma))
        pow[i] <- powCalc(contMat, N, delta = delta, cVal = cVal,
                          twoSide = twoSide, control = control)
      }
      powMat[j,] <- pow
      j <- j + 1
    }
    nams <- dimnames(muMat)[[2]]

    if(!is.character(sumFct)){
      if(length(sumFct) > 1){
        stop("sumFct needs to be a character vector")
      } else {
        sumFct <- deparse(substitute(sumFct)) # for back-compatibility
      }
    }
    nSfunc <- length(sumFct)
    for(i in 1:nSfunc){
      powMat <- cbind(powMat, apply(powMat[,1:nmod, drop = FALSE], 1, sumFct[i] ,...))
    }

    dimnames(powMat) <- list(nSeq, c(nams, sumFct))
    attr(powMat, "alRatio") <- alRatio
    attr(powMat, "sumFct") <- sumFct
    oldClass(powMat) <- "powerMM"
    powMat
  }

### doses <- c(0,10,25,50,100,150)
### models <- list(linear = NULL, emax = 25,
###                logistic = c(50, 10.88111), exponential= 85,
###                betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=T, nrow=2))
### pM <- powerMM(models, doses, base = 0, maxEff = 0.4, sigma = 1, alpha = 0.05,
###                lower = 10, upper = 100, step = 10, scal = 200)
### plot(pM, models="none", line.at = 0.8)

plot.powerMM <-
  function(x, superpose = TRUE, line.at = NULL, models = "all",
           summ = NULL, perc = FALSE, xlab = NULL,
           ylab = ifelse(perc, "Power (%)", "Power"), ...)
  {
    ## Trellis plot of power matrix produced by powerMM
    ## x - matrix with power values for different sample sizes and models
    ## superpose - logical variable indicating if lines should be superpose
    ## line.at - a value, or a vector of values, between 0 and 1, to be
    ##           drawn as horizontal line on plot; maybe most interesting if it
    ##           is set to the desired power value(s) (default: not drawn)
    ## models - which of available models should be included in plot
    ##          "all" and "none" are accepted, else names (or numbers) of models
    ## summ  - summaries to be included in plot; "summary" can be used for
    ##         summary function in x - min and max are also accepted (and
    ##         in combination)
    ## perc - logical indicating if power values should be in percentage
    ## xlab, ylab - labels for x- and y-axis

    ## checking for unbalanced sample sizes
    unbN <- (length(unique(attr(x, "alRatio"))) > 1)
    if (is.null(xlab)) {
      if (unbN) {
        if(sum(attr(x, "alRatio"))==1)
          xlab <- "Overall sample size (unbalanced allocation)"
        else
          xlab <- "Sample size for smallest arm (unbalanced)"
      } else {
        xlab <- "Sample size per dose (balanced)"
      }
    }
    if (perc) x <- 100 * x
    nSeq <- as.integer(dimnames(x)[[1]])
    nams <- dimnames(x)[[2]]
    ## separating model data from summary data
    sumFct <- attr(x, "sumFct")      # original summary functions
    ## default is to use same summary functions as in x
    if (is.null(summ)) summ <- sumFct
    if (!is.null(sumFct)) {
      ## summary functions available in data and requested in summ
      indS <- !is.na(match(sumFct, summ))
      if (any(indS)) {
        summData <- x[, sumFct[indS], drop = FALSE]
      } else {
        summData <- NULL
      }
      modData <- x[, is.na(match(nams, sumFct)), drop = FALSE]
    } else {
      summData <- NULL
      modData <- x
    }
    namsMod <- dimnames(modData)[[2]]
    ## additional summary functions possibly requested in summ
    indSX <- is.na(match(summ, sumFct))
    if(any(indSX)) {
      summ <- summ[indSX]
      if((length(summ) == 1) && (summ == "none")) {
        summDataX <- NULL
      } else {
        nSX <- length(summ)
        summDataX <- vector("list", nSX)
        names(summDataX) <- summ
        for(i in 1:nSX) {
          ## should allow additional arguments to summ[i], but too messy already
          summDataX[[i]] <- apply(modData, 1, summ[i])
        }
        summDataX <- data.frame(summDataX)
      }
    } else {
      summDataX <- NULL
    }
    ## combining the two summary datasets
    if(length(summData) > 0){
      summData <- cbind(summData, summDataX)
    } else {
      summData <- summDataX
    }
    ## checking subset of models
    if (length(models) == 1) {
      if (models != "all") {
        if (models == "none") {
          modData <- NULL
        } else {
          if (is.numeric(models)) {
            if (!is.element(models, 1:ncol(modData))) {
              stop("invalid model selected")
            }
          } else {
            if (!is.element(models, namsMod)) {
              stop("invalid model selected")
            }
          }
          modData <- modData[,models, drop = FALSE]
        }
      }
    } else {                              # need to be subset of models
      if (is.numeric(models)) {
        if (!all(is.element(models, 1:ncol(modData)))) {
          stop("invalid models selected")
        }
      } else {
        if (!all(is.element(models, namsMod))) {
          stop("invalid model selected")
        }
      }
      modData <- modData[,models, drop = FALSE]
    }
    if (length(modData) == 0) x <- summData
    else if (length(summData) == 0) x <- modData
    else x <- cbind(modData, summData)
    x <- data.frame(x)
    nams <- names(x)
    nC <- ncol(x)
    pMatTr <-
      data.frame(pow = as.vector(unlist(x)),
                 n = rep(nSeq, nC),
                 type = factor(rep(nams, each = length(nSeq)), levels = nams))
    if (superpose) {
      panelFunc <- if (is.null(line.at)) {
        panel.superpose
      } else {
        function(x, y, subscripts, groups, lineAt, ...) {
          panel.superpose(x, y, subscripts, groups, ...)
          panel.abline(h = lineAt, lty = 2)
        }
      }
      trLn <- trellis.par.get("superpose.line")[c("col", "lwd", "lty")]
      for(i in seq(along = trLn)) {
        if(length(trLn[[i]]) > nC) trLn[[i]] <- trLn[[i]][1:nC]
      }
      xyplot(pow ~ n, pMatTr, groups = pMatTr$type, subscripts = TRUE,
             panel = panelFunc, type = "l", lineAt = line.at,
             xlab = xlab, ylab = ylab,
             key = list(lines = trLn, text = list(lab = nams), transparent = TRUE,
                        columns = ifelse(nC < 5, nC, min(4,ceiling(nC/min(ceiling(nC/4),3))))), ...)
    } else {                              # models in different panels
      if (is.null(line.at)) {
        panelFunc <- panel.xyplot
      }
      else {
        panelFunc <-
          function(x, y, lineAt, ...) {
            panel.xyplot(x, y, ...)
            panel.abline(h = lineAt, lty = 2) ## used 2 for consistency with above
          }
      }
      xyplot(pow ~ n | type, pMatTr, panel = panelFunc,
             type = "l", lineAt = line.at,
             xlab = xlab, ylab = ylab,
             strip = function(...) strip.default(..., style = 1), ...)
    }
  }

### example call
### doses <- c(0,10,25,50,100,150)
### models <- list(linear = NULL, emax = c(25),
###                logistic = c(50, 10.88111), exponential=c(85),
###                betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=T, nrow=2))
### pM <- powerMM(models, doses, base = 0, maxEff = 0.4, sigma = 1,
###                lower = 10, upper = 100, step = 5, scal = 200)
### plot(pM)
### # reproduces plot in JBS 16, p.651
### plot(pM, line.at = 0.8, models = "none")

sampSize <-
  function(models, doses, base, maxEff, sigma, upperN, lowerN = floor(upperN/2),
           power = 0.8, alRatio = NULL, sumFct = mean, off = 0.1*max(doses), scal = 1.2*max(doses),
           alpha = 0.025, twoSide = FALSE, tol = 1e-3, verbose = FALSE,
           control = mvtnorm.control(), muMat = NULL, typeN = c("arm", "total"), ...)
  {
    ## Calculates  smallest sample size necessary to achieve given power
    ## using a bisection search
    ## models - either a list with model parameters as in planMM
    ##          or an object of class fullMod
    ## doses, maxEff, base, off, scal - as in getPars function
    ## sigma - expected standard deviation
    ## lowerN, upperN - lower and upper limit for the necessary sample size
    ## power - desired power level (between 0 and 1)
    ## alRatio - allocation ratios ()
    ## sumFct - summary function for the individual power values
    ## alpha - level of significance
    ## twoSide - logical variable indicating if one or two-sided test
    ## tol - tolerance for power calculation
    ## verbose - logical controlling amount of output in iterations
    ## control - see mvtnorm.control function
    ## typeN - type of sample size specified in lowerN and upperN, in the case of
    ##         unequal allocation: "arm" refers to smallest trt arm, "total" to
    ##         overall sample size
    ## ... - possible additional arguments for sumFct

    if (missing(models)) {
      ## may pass necessary information via muMat
      if (is.null(muMat)) {
        stop("either models or muMat need to be specified")
      }
    } else {
      ## if models is of class fullMod, will ignore other arguments
      if (inherits(models, "fullMod")) {
        doses <- attr(models, "doses")
        base <- attr(models, "base")
        maxEff <- attr(models, "maxEff")
        off <- attr(models, "off")
        scal <- attr(models, "scal")
        Models <- models
      } else {
        Models <- fullMod(models, doses, base, maxEff, off, scal)
      }
      muMat <- modelMeans(Models, doses, FALSE, off, scal)
    }
    nD <- length(doses)
    Ntype <- match.arg(typeN)
    if(!is.null(alRatio)){
      if(length(alRatio)!= nrow(muMat)){
        stop("alRatio needs to be of same length as number of dose levels")
      }
      if(any(alRatio <= 0)){
        stop("all entries of alRatio need to be positive")
      }
      if (Ntype == "arm") {
        alRatio <- alRatio/min(alRatio)
      } else {
        alRatio <- alRatio/sum(alRatio)
      }
    } else {
      alRatio <- 1
    }
    contMat <- modContr(muMat, upperN * alRatio)       # no need to change
    summPow <-
      function(n, nD, contMat, muMat, sigma, alpha, sumFct, control, twoSide, ...)
      {
        ## returns "combined" power for a specified sample size n
        pow <- numeric(ncol(muMat))
        cVal <- critVal(contMat, n, alpha, control, twoSide)
        for(i in 1:ncol(muMat)) {
          pow[i] <- powCalc(contMat, n, cVal = cVal, mu = muMat[,i],
                            sigma = sigma, twoSide = twoSide,
                            control = control)
        }
        res <- sumFct(pow, ...)
        names(pow) <- dimnames(contMat)[[2]]
        attr(res, "power under models") <- round(pow, 4)
        res
      }
    ## need to ensure that limits on n bracket power
    upperPow <-
      summPow(round(upperN * alRatio), nD, contMat, muMat,
              sigma, alpha, sumFct, control, twoSide, ...)
    if(upperPow < power){
      if(Ntype == "arm") sent <- "sample size in smallest group"
      else sent <- "total sample size"
      warning("upper limit for ", sent," is raised")
    }
    while (upperPow < power) {
      upperN <- 2 * upperN
      upperPow <-
        summPow(round(upperN * alRatio), nD, contMat, muMat,
                sigma, alpha, sumFct, control, twoSide, ...)
    }
    lowerPow <-
      summPow(round(lowerN * alRatio), nD, contMat, muMat, sigma,
              alpha, sumFct, control, twoSide, ...)
    if(lowerPow > power){
      if(Ntype == "arm") sent <- "sample size in smallest group"
      else sent <- "total sample size"
      warning("lower limit for ", sent," is decreased")
    }
    while (lowerPow > power) {
      lowerN <- round(lowerN/2)
      if (lowerN == 0) stop("cannot find lower limit on n")
      lowerPow <-
        summPow(round(lowerN * alRatio), nD, contMat, muMat, sigma,
                alpha, sumFct, control, twoSide, ...)
    }
    ## now use simple binary search
    if (verbose) {
      cat("Upper N:", upperN, "Upper power", round(upperPow, 3), "\n")
      cat("Lower N:", lowerN, "Lower power", round(lowerPow, 3), "\n\n")
    }
    currPow <- 0
    niter <- 0
    while (abs(currPow - power) > tol & (upperN > lowerN + 1)) {
      currN <- round((upperN + lowerN)/2)
      currPow <-
        summPow(round(currN * alRatio), nD, contMat, muMat, sigma,
                alpha, sumFct, control, twoSide, ...)
      if (currPow > power) {
        upperN <- currN
      } else {
        lowerN <- currN
      }
      niter <- niter + 1
      if (verbose){
        cat("Iter: ", niter, ", N = ", currN, ", power = ", round(currPow, 3),
            "\n", sep = "")
      }
    }
    while (currPow < power) { # step up until currPow >= power
      currN <- currN + 1
      currPow <-
        summPow(round(currN * alRatio), nD, contMat, muMat, sigma,
                alpha, sumFct, control, twoSide, ...)
    }
    if(length(alRatio) > 1) {   # unequal allocation
      if (Ntype == "arm") {
        res <- list(samp.size = round(currN * alRatio),
                    approx.power = round(currPow, 4))
      } else {
        res <- list(samp.size = round(currN * alRatio),
                    approx.power = round(currPow, 4))
      }
    } else {
      res <- list(samp.size = currN, approx.power = round(currPow, 4))
    }
    attr(res, "alRatio") <- round(alRatio/min(alRatio), 4)
    attr(res, "power") <- power
    attr(res, "twoSide") <- twoSide
    attr(res, "alpha") <- alpha
    attr(res, "sumFct") <- deparse(substitute(sumFct))
    oldClass(res) <- "sampSize"
    res
  }

print.sampSize <- function(x, ...){
  cat("MCPMod sampSize\n")
  cat("\n")
  cat("Input parameters:\n")
  cat(" Summary Function:", attr(x, "sumFct"), "\n")
  cat(" Desired combined power value:", attr(x, "power"),"\n")
  if(attr(x, "twoSide")) side <- "two-sided"
  else side <- "one-sided"
  cat(" Level of significance:", paste(attr(x, "alpha"), " (", side, ")", sep=""), "\n")
  alRatio <- attr(x, "alRatio")
  if(length(alRatio) == 1){
    alRatio <- "balanced"
    sampMsg <- "Sample size per group:"
  } else {
    sampMsg <- "Sample sizes:"
  }
  cat(" Allocations:", alRatio, "\n\n")
  cat(paste(sampMsg), paste(x$samp.size), "\n")
  cat("\n")
  cat("Associated", attr(x, "sumFct"),"power:", paste(x$approx.power),"\n")
  cat("Power under models:","\n")
  print(attr(x$approx.power, "power under models"))
  cat("\n")
}

### example call
### doses <- c(0,10,25,50,100,150)
### models <- list(linear = NULL, emax = c(25),
###                logistic = c(50, 10.88111), exponential=c(85),
###                betaMod=matrix(c(0.33,2.31,1.39,1.39), byrow=T, nrow=2))
### sampSize(models, doses, base = 0, maxEff = 0.4, sigma = 1,
###          alpha = 0.05, upperN = 100, scal = 200)
### sampSize(models, doses, base = 0, maxEff = 0.4, sigma = 1,
###        alpha = 0.05, upperN = 400, scal = 200, alRatio=c(2,1,1,1,1,1))
### dvec <- c(0, 10, 50, 100)
### mu1 <- c(1, 2, 2, 2)
### mu2 <- c(1, 1, 2, 2)
### mu3 <- c(1, 1, 1, 2)
### mMat <- cbind(mu1, mu2, mu3)
### dimnames(mMat)[[1]] <- dvec
### sampSize(muMat = mMat, doses = dvec, sigma = 1,
###        alpha = 0.05, upperN = 400, alRatio=c(2,1,1,1))


LP <- function(models, model, type = c("both", "LP1", "LP2"), paramRange, doses,
               base, maxEff, sigma, n,
               len = c(10, 1), nr = 1, alpha = 0.025, twoSide = FALSE,
               off = 0.1*max(doses), scal = 1.2*max(doses), control = mvtnorm.control()){
  ## function to calculate loss in power associated with
  ## prior parameter mis-specification
  ## models, doses, base, maxEff, sigma as above
  ## model - model for which LP should be investigated
  ## type - one of LP1 or LP2 (see JBS paper for definition)
  ## paramRange - vector of lower and upper limit for the prior parameter
  ##       for models with two prior parameters matrix with the ranges
  ##       for each prior parameter in the rows
  ## n - number (in case of balanced sample size)
  ##      or vector of group sample sizes
  ## len - number of points between min(paramRange) and max(paramRange)
  ##       on which LP is calculated has to be of length 2 in case of
  ##       models with 2 prior parameter
  ## nr - in case there is more than one model from one model class
  ##   in the candidate set (e.g. two emax models), nr gives the number of the model
  ## alpha, twoSide, scal, off as above

  type <- match.arg(type)
  builtIn <- c("emax", "exponential", "quadratic",
               "betaMod", "logistic", "sigEmax")
  usedM <- names(models)
  allowedMod <- intersect(builtIn, usedM)
  if(!is.element(model, allowedMod)){
    msg <- paste("model must be one of", paste(allowedMod, collapse=", "), "\n")
    stop(msg)
  }
  if (inherits(models, "fullMod")) {
    doses <- attr(models, "doses")
    base <- attr(models, "base")
    maxEff <- attr(models, "maxEff")
    off <- attr(models, "off")
    scal <- attr(models, "scal")
    Models <- models
  } else {
    Models <- fullMod(models, doses, base, maxEff, off, scal)
  }
  nD <- length(doses)
  if(length(n)==1) n <- rep(n,nD)
  nDF <- sum(n) - nD
  # calculate contrasts and critical value
  pM <- planMM(Models, doses, n, off, scal, FALSE, alpha, twoSide, control)
  contMat <- pM$contMat
  critV <- pM$critVal
  corMat <- pM$corMat
  # calculate mean vector under 'model'
  # determine, whether there is more than one prior estimate for model "model"
  if(nrInd <- length(Models[[model]]) > 4) mod <- list(Models[[model]][nr,])
  else mod <- list(Models[[model]])
  names(mod) <- model
  mu <- modelMeans(mod, doses, std=FALSE, off, scal)
  # adjust arguments
  twoPars <- is.element(model, c("logistic", "betaMod", "sigEmax"))
  if(!twoPars){ # put paramRange in the same format
    paramRange <- matrix(c(paramRange, 0, 0), nrow=2, byrow=TRUE)
    len <- c(len,1)
  }
  if(type[1] != "LP1"){ # "LP2" or "both"
    ModelsIt <- Models
  }
  if(type[1] != "LP2"){ # "LP1" or "both"
    # calculate nominal power
    nomPower <- powCalc(contMat, n, alpha, mu = mu, sigma = sigma,
                        cVal = critV, corMat = corMat,
                        twoSide = twoSide, control = control)[1]
  }
  actPower <- potPower <- numeric(prod(len))
  j <- 1
  # parameter ranges
  sq1 <- seq(min(paramRange[1,]),max(paramRange[1,]),length=len[1])
  sq2 <- seq(min(paramRange[2,]),max(paramRange[2,]),length=len[2])

  for(h in sq2) {
    for(i in sq1) {
      if(twoPars) mod <- list(c(i,h))
      else mod <- list(i)
      names(mod) <- model
      mod <- fullMod(mod, doses, base, maxEff, off, scal)
      # mean vector under model with prior parameter i,h
      mu <- modelMeans(mod, doses, std=FALSE, off, scal)
      actPower[j] <- powCalc(contMat, n, alpha, mu = mu, sigma = sigma,
                             cVal = critV, corMat = corMat,
                             twoSide = twoSide, control = control)
      if(type[1] == "LP2" | type[1] == "both"){  # recalculate opt. cont. and critical value
        if(twoPars){          # using i,h as prior parameters
          if(nrInd) ModelsIt[[model]][nr,3:4] <- c(i,h)
          else ModelsIt[[model]][3:4] <- c(i,h)
        } else {
          if(nrInd) ModelsIt[[model]][nr,3] <- i
          else ModelsIt[[model]][3] <- i
        }
        pMIt <- planMM(ModelsIt, doses, n, off, scal, FALSE, alpha, twoSide)
        potPower[j] <- powCalc(pMIt$contMat, n, alpha, mu = mu, sigma = sigma,
                               cVal = pMIt$critVal, corMat = pMIt$corMat,
                               twoSide = twoSide, control = control)
      }
      j <- j + 1
    }
  }
  # determine parametername for output
  par <- switch(model,
                logistic = c("ED50", "delta"),
                betaMod = c("delta1", "delta2"),
                sigEmax = c("ED50", "h"),
                emax = "ED50",
                "delta")
  p1 <- rep(sq1,len[2])
  if(twoPars){
    p2 <- rep(rep(round(sq2,2),length=len[2]), each=len[1]) # values of 2nd parameter
    if(nrInd) used <- Models[[model]][nr,3:4] # determine used value
    else  used <- Models[[model]][3:4]
  } else {
    if(nrInd) used <- Models[[model]][nr,3]
    else  used <- Models[[model]][3]
  }
  res <- cbind(actPower) # matrix for output
  if(type[1] == "LP1"){
    LP1 <- nomPower - actPower
    res <- cbind(res, LP1)
  } else if(type[1] == "LP2"){
    res <- cbind(potPower, res)
    LP2 <- potPower - actPower
    res <- cbind(res, LP2)
  } else if(type[1] == "both"){
    res <- cbind(potPower, res)
    LP1 <- nomPower - actPower
    LP2 <- potPower - actPower
    res <- cbind(res, LP1, LP2)
  }
  if(twoPars){  # names for prior parameters
    res <- cbind(p1, p2, res)
    dimnames(res)[[2]][1:2] <- par
  } else {
    res <- cbind(p1, res)
    dimnames(res)[[2]][1] <- par
  }
  attr(res, "model") <- model
  attr(res, "len") <- len
  attr(res, "used") <- used
  attr(res, "sampSizes") <- n
  attr(res, "type") <- type[1]
  attr(res, "twoPars") <- twoPars
  if(type[1] != "LP2") attr(res, "nomPower") <- nomPower
  oldClass(res) <- "LP"
  res
}

print.LP <- function(x, digits = 3, ...){
  cat("MCPMod LP","\n")
  cat("\n")
  cat("Model:", paste(attr(x, "model")),"\n")
  cat("Used Prior Parameter:", paste(attr(x, "used")), "\n")
  cat("Sample Sizes:", paste(attr(x, "sampSizes")),"\n")
  type <- attr(x, "type")
  if(type=="LP1" | type=="both"){
    cat("Nominal Power:", paste(round(attr(x, "nomPower"), digits)),"\n")
  }
  cat("\n")
  cat("Calculated Power Values:", "\n")
  print(round(data.frame(unclass(x)), digits = digits))
}

plot.LP <- function(x, line = TRUE, type = NULL, spldf = 5, ...){
  ## line - logical determining whether smooth line
  ##        should be fit to power values
  ## type - one of "LP1","LP2","both" or NULL
  ##        if NULL the type of the LP call is used

  if(line){
    val <- 1
    len <- attr(x, "len")
    if(len[1] < 4)
      stop("at least 4 points needed to obtain smooth curve. Try: line = FALSE")
    if (!(spldf > 1 && spldf <= dim(x)[1]))
      warning("you must supply 1 < spldf <= len")
  }
  twoPars <- attr(x, "twoPars")
  model <- attr(x, "model")
  used <- attr(x, "used")
  par <- dimnames(x)[[2]][1:(1+twoPars)]
  par1 <- x[,1]
  if(is.na(par[2])){
    par2 <- rep(1, dim(x)[1])
  } else {
    par2 <- x[,2]
  }
  if(is.null(type)){  # check if type is given in LP matrix
    type <- attr(x, "type")
  } else {
    if(attr(x, "type")!="both" & attr(x, "type")!=type){
      stop("invalid type selected")
    }
  }
  # matrix for plotting data
  if(type == "both"){
    type <- c("LP1", "LP2")
    grp <- rep(type, each = dim(x)[1])
    pData <- data.frame(LP = c(x[, type]), par1 = rep(par1, 2),
                        par2 = rep(paste(paste(par[2],"="), par2), 2),
                        group = grp)
  } else {
    grp <- rep(type, dim(x)[1])
    pData <- data.frame(LP = x[, type], par1 = par1,
                        par2 = paste(paste(par[2],"="), par2),
                        group = grp)
  }
  main <- list(paste("Model:", model, ", Used value:", paste(used, collapse=" ")))
  L <- length(unique(grp))
  spL <- trellis.par.get("superpose.symbol")
  col <- spL$col[1:L] # use the dot color for the lines
  if(L == 2){ # key only needed when both LP1 and LP2 are to be plotted
    keyList <- list(points = list(col = col, pch = 1), transparent = TRUE,
                    text = list(paste(type), cex = 0.9), columns =  2)
    ylab <- NULL
  } else {
    keyList <- NULL
    ylab <- type
  }
  xyplot(LP~par1|par2, data=pData,
         panel.data = list(uv = used[1], line = line, col = col, df = spldf),
         groups = grp,
         panel = function(x, y, subscripts, groups,..., panel.data) {
           panel.abline(h = 0, lty = 2)
           panel.abline(v = panel.data[["uv"]], lty = 2)
           s <- seq(min(x), max(x), length = 101)
           panel.superpose(x, y, type = "p", lwd = 1.3, groups = groups, subscripts = subscripts)

           j <- 1
           if(panel.data[["line"]]){
             for(i in unique(groups[subscripts])){   # add smooth line
               ind <- which(groups[subscripts]==i)
               p <- predict(smooth.spline(x[ind], y[ind], df = panel.data[["df"]]), s)
               panel.xyplot(p$x, p$y , type = "l", col = panel.data[["col"]][j])
               j <- j + 1
             }
           }
         },
         main = main, xlab=par[1], strip=twoPars, ylab = ylab,
         key = keyList)
}



#######################################################################
## Copyright 2008, Novartis Pharma AG
##
## This program is Open Source Software: you can redistribute it
## and/or modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation, either version 3 of
## the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program. If not, see http://www.gnu.org/licenses/.

## functions to analyze dose finding data according to MCP-Mod

## t-stats for MCP step
getTstat <-
  function(data, n, contMat, resp = "resp", dose = "dose")
  {
    if (any(is.na(match(c(resp, dose), names(data))))) {
      stop(resp," and/or ", dose, " not found in data")
    }
    dose <- data[, dose]
    resp <- data[, resp]
    ## means per dose
    mn <- tapply(resp, dose, mean)

    ## pooled standard deviation
    sdv <- tapply(resp, dose, sd)
    if (length(n) == 1) n <- rep(n, length(sdv))
    # remove NAs for dose groups with only 1 obs.
    if(any(n==1)){
      sdv[n==1] <- 0
    }
    n1 <- n - 1
    sdv <- sqrt(sum(n1 * sdv^2)/sum(n1))

    ## contrasts
    ct <- as.vector(mn %*% contMat)
    den <- sdv * sqrt(colSums((contMat^2)/n))

    ## max t-stat
    val <- ct/den
    names(val) <- dimnames(contMat)[[2]]
    val
  }

# function to calculate p-values
pValues <-
  function(cMat, n, alpha = 0.025, Tstats, control = mvtnorm.control(),
           twoSide = FALSE, corMat = NULL, ...)
  {
    ## function to calculate p-values
    ##
    ## cMat - model contrast matrix nDose x nMod
    ## n - scalar or vector with sample size(s)
    ## alpha - significance level
    ## control - see mvtnorm.control function
    ## twoSide - logical variable indicating if one or two-sided test
    ## corMat - correlation matrix

    nD <- nrow(cMat)
    nMod <- ncol(cMat)
    if(length(Tstats) != nMod){
      stop("Tstats needs to have length equal to the number of models")
    }
    if (length(n) == 1) {
      nDF <- nD * (n - 1)
    } else {
      if (length(n) != nD) {
        stop("'n' must have length as number of doses")
      }
      nDF <- sum(n) - nD
    }
    if(is.null(corMat)){
      corMat <- t(cMat)%*%(cMat/n)
      den  <- sqrt(crossprod(t(colSums(cMat^2/n))))
      corMat <- corMat / den
    }
    ctrl <- mvtnorm.control()
    if (!missing(control)) {
      control <- as.list(control)
      ctrl[names(control)] <- control
    }
    if(twoSide){
      lower <- matrix(rep(-Tstats, each = nMod), nrow = nMod)
    }
    else lower <- matrix(rep(-Inf, nMod^2), nrow = nMod)
    upper <- matrix(rep(Tstats, each = nMod), nrow = nMod)
    pVals <- numeric(nMod)
    for(i in 1:nMod){
      pVals[i] <-   1 - pmvt(lower[,i], upper[,i], df = nDF, corr = corMat,
                             algorithm = ctrl, ...)
    }
    pVals
  }

### functions to fit models, and functions to get starting values
### for the nls function

### Initial estimates for built-in models
### two auxiliary functions to get starting values for the nls algorithm
getInitLRasymp <-
  ## left and right asymptotes
  function(meanVal, dlt = 0.1)
  {
    rg <- range(meanVal)
    dlt <- dlt * diff(rg)
    c(e0 = rg[1] - dlt, eMax = rg[2] + dlt)
  }

getInitP <-
  ## dose that gives certain percentage of max effect
  function(meanVal, eVec = getInitLRasymp(meanVal), p = 0.5, dose, ve=FALSE)
  {
    e0   <- eVec[1]
    eMax <- eVec[2]
    targ <- e0 * (1 - p) + eMax * p
    ind  <- meanVal <= targ
    ind1 <- !ind
    if (!any(ind) || all(ind)) stop("cannot find initial values")
    ind  <- ind & (dose <= max(dose[ind1]))
    d1   <- max(dose[ind]); p1 <- meanVal[dose == d1]
    ind1 <- ind1 & dose > d1
    d2   <- min(dose[ind1]); p2 <- meanVal[dose == d2]
    res <- d1 + (d2 - d1) * (targ - p1)/(p2 - p1)
    names(res) <- NULL
    res
  }

getInitBeta <-
  function(m, scal, dose)
  {
    ## calculates starting values for the beta model
    ## for the nls algorithm
    ## dose and m assumed to be ordered
    ## in the same order (according to dose)
    dmax <- dose[which(m==max(m))]
    e0 <- m[dose==0]
    emax <- max(m) - e0
    ds <- dose[order(m)[length(m)-1]] # select as 2nd dose the one
    # with 2nd highest resp.
    z <- dmax/(scal-dmax)
    f <- (m[dose==ds] - e0)/emax
    beta <- try(log(f)/(log(ds^z*(scal-ds)) - log(dmax^z*(scal-dmax))))
    alpha <- z*beta
    if(is.na(alpha)) alpha <- beta <- 1 # may occur if f < 0; (1,1) as
    # initial values seems always to
    # work quite well
    res <- c(alpha, beta)
    names(res) <- NULL
    res
  }

getInit <-
  ## Initial estimates for built-in models, if needed
  ## non-linear fitting is done with p-linear, but Gauss-Newton used
  ## later due to some problems with se.fit of the predictions of the
  ## p-linear fit
  function(data, model = c("emax", "exponential", "logistic", "betaMod",
                           "sigEmax"), scal)
  {
    model <- match.arg(model)             # only built-in allowed
    meanVal <- data$respM
    dose <- data$doseM

    switch(model,
           emax = {
             c(ed50 = getInitP(meanVal, dose = dose))
           },
           exponential = {
             e0 <- getInitLRasymp(meanVal)[1]
             meanVal <- meanVal - e0
             aux <- coef(lm(log(meanVal) ~ dose, na.action = na.omit))
             names(aux) <- NULL
             c(delta = 1/aux[2])
           },
           logistic = {
             ed50 <- getInitP(meanVal, dose = dose)
             delta <- getInitP(meanVal, p = 0.75, dose = dose) - ed50
             c(ed50 = ed50, delta = delta)
           },
           betaMod = {
             init <- getInitBeta(meanVal, scal, dose)
             c(delta1 = init[1], delta2 = init[2])
           },
           sigEmax = {
             ed50 <- getInitP(meanVal, dose = dose)
             ed10 <- getInitP(meanVal, p = 0.1, dose = dose)
             ed90 <- getInitP(meanVal, p = 0.9, dose = dose)
             h <- 1.91/log10(ed90/ed10) # From Ch. 14, Ting (2006),
             # Dose finding in drug development
             c(ed50 = ed50, h = h)
           })
  }

fitModel <-
  ## fits dose-response model from list of built-in models,
  ## or according to model provided by user
  function(data, model, resp = "resp", dose = "dose", start,
           control = NULL, off = 1, scal = 1, uModPars = NULL,
           addArgs = NULL, na.action = na.fail)
  {
    ## data - data.frame with information on dose and response
    ##          and possibly other covariates used in model
    ## model - string with name of model function to be used
    ## resp, dose - characters with names of columns in data
    ##              corresponding to response and dose
    ## control - optional list with control parameters for fit
    ## off - offset for dealing with dose = 0 in lin-log model
    ## uModPars - optional character vector with names/expressions
    ##            of user-defined model parameters (names(start) used by
    ##            default)
    ## addArgs - optional character vector with names of additional
    ##           arguments (variables) to user-defined model
    ## na.action - function specifying how to handle NAs in data

    ## ensuring resp and dose are defined in data
    data$dose <- data[, dose]
    data$resp <- data[, resp]

    ## reacting to possible NAs
    ## it would be better to only take into account NAs in
    ## varibles actually used in the model.
    data <- na.action(data)

    ## checking if built-in model and assigning it number, if so
    builtIn <- c("linlog", "linear", "quadratic", "emax",
                 "exponential","logistic", "betaMod", "sigEmax")
    modelNum <- match(model, builtIn)

    ## under the built-in models with no covariates, the mean
    ## modeling can be used in all cases
    respM <- tapply(data$resp, data$dose, mean) # fit models based on means
    doseM <- sort(unique(data$dose))
    n <- as.vector(table(data$dose))
    dataM <- data.frame(doseM = doseM, respM = respM, n = n)
    vars <- tapply(data$resp, data$dose, var)
    # allow for dose groups with only one patient
    vars[n==1] <- 0 # replace NAs with 0
    S2 <- sum((n - 1) * vars)

    if (!is.na(modelNum)) {               # built-in model
      if (is.element(modelNum, 1:3)) {    # linear models
        if (modelNum == 1) {
          dataM$off <- rep(off, nrow(dataM))
          ## linear-log model
          fm <- lm(respM ~ I(log(doseM + off)), dataM, qr=TRUE, weights = n)
          base <- coef(fm)[1]+coef(fm)[2]*log(off)
        }
        if (modelNum == 2) {
          ## linear model
          fm <- lm(respM ~ doseM, dataM, qr=TRUE, weights = n)
          base <- coef(fm)[1]
        }
        if (modelNum == 3) {
          ## quadratic model
          fm <- lm(respM ~ doseM + I(doseM^2), dataM, qr=TRUE, weights = n)
          base <- coef(fm)[1]               # baseline
        }
      } else {                            # nonlinear built-in model
        if (is.null(start)) {             # need to derive initial values
          start <- getInit(dataM, model, scal)
        }
        if (modelNum == 4) {      # Emax model
          names(start) <- NULL
          start <- c(led50 = log(start))
          fm <- try(nls(respM ~ cbind(1, emax(doseM, 0, 1, exp(led50))),
                        start = start, data = dataM, model = TRUE,
                        algorithm = "plinear", control = control,
                        weights = n), silent = TRUE)
          if (!inherits(fm, "nls")) {
            fm <- NA
          } else {
            base <- coef(fm)[2]
          }
        }
        if (modelNum == 5) {               # exponential model
          names(start) <- NULL
          start <- c(ldelta = log(start))
          fm <-
            try(nls(respM ~ cbind(1, exponential(doseM, 0, 1, exp(ldelta))),
                    start = start, weights = n, data = dataM, model = TRUE,
                    control = control, algorithm = "plinear"), silent = TRUE)
          if (!inherits(fm, "nls")) {
            fm <- NA
          } else {
            base <- coef(fm)[2]
          }
        }
        if (modelNum == 6) {               # logistic model
          start <- c(log(start["ed50"]), start["delta"])
          names(start) <- c("led50", "delta")
          fm <- try(nls(respM ~ cbind(1, logistic(doseM, 0, 1, exp(led50),
                                                  delta)),
                        start = start, algorithm = "plinear", data = dataM,
                        model = TRUE, control = control, weights = n),
                    silent = TRUE)
          if (!inherits(fm, "nls")) {
            fm <- NA
          } else {
            base <- predict(fm, newdata = list(doseM = 0))
          }
        }
        if (modelNum == 7) {              # beta model
          dataM$scal <- rep(scal, nrow(dataM))
          fm <-
            try(nls(respM ~ cbind(1, betaMod(doseM, 0, 1, delta1, delta2, scal)),
                    start = start, weights = n, data = dataM, model = TRUE,
                    algorithm = "plinear", control = control), silent = TRUE)
          if (!inherits(fm, "nls")) {
            fm <- NA
          } else {
            base <- coef(fm)[3]
          }
        }
        if (modelNum == 8) {              # sigEmax model
          start <- c(log(start["ed50"]), start["h"])
          names(start) <- c("led50", "h")
          fm <-
            try(nls(respM ~ cbind(1, sigEmax(doseM, 0, 1, exp(led50), h)),
                    start = start, model = TRUE,
                    algorithm = "plinear", control = control, weights = n),
                silent = TRUE)
          if (!inherits(fm, "nls")) {
            fm <- NA
          } else {
            base <- coef(fm)[3]
          }
        }
      }
      if (length(fm) > 1) ResidSS <- S2 + deviance(fm)
    } else {                              # user defined model
      ## first check for starting estimates, parameter names, etc
      if (is.null(start)) {
        stop("must provide stating estimates for user-defined model")
      }
      namStart <- names(start)
      if (is.null(namStart)) {
        stop("'start' must have names for user-defined models")
      }
      if (is.null(uModPars)) uModPars <- namStart  # default parameters
      ## create the model call
      modForm <- paste("resp ~ ", model, "(dose,",
                       paste(uModPars, collapse = ","), sep = "")
      if (!is.null(addArgs)) {
        modForm <- paste(modForm, ",", paste(addArgs, collapse = ","))
      }
      modForm <- paste(modForm, ")")
      modForm <- eval(parse(text = modForm))
      ## will assume non-linear model
      fm <- try(do.call("nls", list(modForm, data, start, control)))
      if (!inherits(fm, "nls")) {
        fm <- NA
      } else {
        ResidSS <- deviance(fm)
        dataM <- NULL
        ## following will not work when covariates are used in nls model
        base <- predict(fm, newdata = list(dose = 0))
      }
    }
    if (length(fm) > 1) {
      res <- list(fm = fm, base = base)
      ## saving information for other methods
      attr(res, "ResidSS") <- ResidSS
      attr(res, "model") <- model
      attr(res, "scal") <- scal
      attr(res, "off") <- off
      attr(res, "dataM") <- dataM
    } else {
      res <- NA
    }
    oldClass(res) <- "fitMCPMod"          # allow use of methods later
    res
  }

# function to return analyt. gradient for built-in or
# user defined models
getGrad <-
  function(obj, dose, uGrad = NULL)
    ## obj - an object inheriting from class fitMCPMod
    ## dose - vector with doses at which to calculate the  gradient
    ## uGrad - optionally, a character string with the name of a
    ##         user-defined gradient function
  {
    if(!inherits(obj, "fitMCPMod")) {
      stop("obj must inherit from class fitMCPMod")
    }
    model <- attr(obj, "model")
    fm <- obj$fm
    if (!inherits(fm, "nls")) {
      ## linear models are left out
      stop("fitted object must inherit from class nls")
    }
    ## only used internally in getGrad
    lg2 <- function(x) ifelse(x == 0, 0, log(x))

    cf <- coef(obj)                       # estimated parameters
    res <-
      switch(model,
             "emax" = {
               eMax <- cf[2]; ed50 <- cf[3]
               cbind(1, dose/(ed50 + dose), -eMax*dose/(dose + ed50)^2)
             },
             "logistic" = {
               eMax <- cf[2]; ed50 <- cf[3]; delta <- cf[4]
               den <- 1 + exp((ed50 - dose)/delta)
               g1 <- -eMax*(den-1)/(delta*den^2)
               g2 <- eMax*(den-1)*(ed50-dose)/(delta^2*den^2)
               cbind(1, 1/den, g1, g2)
             },
             "sigEmax" = {
               eMax <- cf[2]; ed50 <- cf[3]; h <- cf[4]
               den <- (ed50^h + dose^h)
               g1 <- dose^h/den
               g2 <- -ed50^(h-1)*dose^h*h*eMax/den^2
               g3 <- eMax*dose^h*ed50^h*lg2(dose/ed50)/den^2
               cbind(1, g1, g2, g3)
             },
             "betaMod" = {
               scal <- attr(obj, "scal")
               dose <- dose/scal
               if (any(dose > 1)) {
                 stop("doses cannot be larger than scal in betaModel")
               }
               delta1 <- cf[3]; delta2 <- cf[4]; eMax <- cf[2]
               maxDens <- (delta1^delta1)*(delta2^delta2)/
                 ((delta1 + delta2)^(delta1+delta2))
               g1 <- ((dose^delta1) * (1 - dose)^delta2)/maxDens
               g2 <- g1*eMax*(lg2(dose)+lg2(delta1+delta2)-lg2(delta1))
               g3 <- g1*eMax*(lg2(1-dose)+lg2(delta1+delta2)-lg2(delta2))
               cbind(1, g1, g2, g3)
             },
             "exponential" = {
               delta <- cf[3]
               e1 <- cf[2]
               cbind(1, exp(dose/delta) - 1, -exp(dose/delta)*dose*e1/delta^2 )
             },
             {
               ## user defined gradient function
               if(is.null(uGrad)) {
                 stop("user-defined gradient needs to be specified")
               }
               do.call(uGrad, c(list(dose), cf))
             })
    res
  }

# function to return coefficients from fit
coef.fitMCPMod <-
  function(object, ...)
  {
    ## object - object inheriting from class MCPMod
    fm <- object$fm
    cf <- coef(fm)
    if(inherits(fm, "lm")) return(cf)
    temp <- names(cf) # save names for user-models
    names(cf) <- NULL
    model <- attr(object, "model")
    cf <- switch(model,
                 "emax" = c(e0 = cf[2], eMax = cf[3], ed50 = exp(cf[1])),
                 "logistic" = c(e0=cf[3], eMax=cf[4], ed50 = exp(cf[1]),
                                delta = cf[2]),
                 "exponential" = c(e0 = cf[2], e1 = cf[3], delta = exp(cf[1])),
                 "betaMod" = c(e0=cf[3], eMax=cf[4], delta1=cf[1], delta2=cf[2]),
                 "sigEmax" = c(e0=cf[3], eMax=cf[4], ed50=exp(cf[1]), h=cf[2]),
                 {
                   names(cf) <- temp
                   cf
                 }
    )
    cf
  }

# predict method
predict.fitMCPMod <-
  function(object, doseSeq, se.fit = TRUE, uGrad = NULL, ...)
  {
    # object - either a lm object or a nls object in the latter case
    #     the code distinguishes between user-defined models
    #     and built-in models
    # doseSeq - sequence for predictions
    # se.fit - logical indicating whether sd of the mean should be
    #          calculated
    model <- attr(object, "model")
    scal <- attr(object, "scal")
    off <- attr(object, "off")
    fm <- object$fm
    if(inherits(fm, "lm")){ # in this case use the implemented predict method
      newdata <- data.frame(doseM = doseSeq, off = rep(off, length(doseSeq)))
      res <- predict(fm, newdata, se.fit = se.fit)
      ## Need to correct se.fit, if requested
      if(se.fit) {
        df <- sum(fm$weights) - length(coef(fm))
        sigCorr <- sqrt(attr(object, "ResidSS")/df)
        sigOrig <- summary(fm)$sigma
        res$se.fit <- (sigCorr/sigOrig) * res$se.fit
        res$df <- df
        res$residual.scale <- sigCorr
      }
      return(res)
    }
    builtIn <- is.element(model, c("emax", "exponential", "logistic",
                                   "betaMod", "sigEmax"))
    if(inherits(fm, "nls") & !builtIn){ # user defined model
      mn <- predict(fm, data.frame(dose = doseSeq))
      if(se.fit) {
        ## first, check if model includes a gradient attribute
        grd <- attr(mn, "gradient")
        if(is.null(grd)) {
          if(is.null(uGrad)) {
            stop("neither gradient attribute, nor uGrad defined")
          }
          grd <- getGrad(object, doseSeq, uGrad)
        }
        Rinv <- solve(fm$m$Rmat())
        sig <- summary(fm)$sigma
        seFit <- sig * sqrt(rowSums((grd %*% Rinv)^2))
        df <- df.residual(fm)
        res <- list(fit = mn, se.fit =  as.vector(seFit),
                    residual.scale = sig, df = df)
        return(res)
      } else {
        return(mn)
      }
    } else { # built in model
      dd <- data.frame(doseM = doseSeq) # dose levels to be predicted
      cf <- coef(object)  # extract coefficients
      if(model == "betaMod"){
        cf <- c(cf, scal) # add scale parameter for betaModel
        dd$scal <- scal
      }
      mn <- predict(fm, dd) # predict mean response
      if(!se.fit){
        return(mn)
      }
      RSS <- attr(object, "ResidSS")   # get residual sum of squares
      doseV <- attr(object, "dataM")$doseM # recover dose-vector
      n <- fm$weights
      df <- sum(n) - length(cf)
      sig <- sqrt(RSS/df) # residual variance estimate
      ## Gradient matrix: see Bates/Watts p. 58/59
      V <- getGrad(object, doseV)
      if(all(!is.infinite(V)) & all(!is.nan(V))){
        ## checking for infinities (can happen because the analytic
        ## gradient is not used for fitting)
        V <- t(crossprod(V, sqrt(diag(n))))  # weights for unequal sample sizes
        R <- qr.R(qr(V))
        Rinv <- try(solve(R))
        if (!inherits(Rinv, "matrix")){
          ## sometimes R is singular (typically for very unequally
          ## spaced dose designs)
          warning("Cannot cannot calculate standard deviation for ",
                  model, " model.\n")
          seFit <- rep(NA, length(doseSeq))
        } else {
          v <- getGrad(object, doseSeq) # calc. gradient
          seFit <- sig * sqrt(rowSums((v%*%Rinv)^2))
        }
      } else {
        warning("Cannot cannot calculate standard deviation for ",
                model, " model.\n")
        seFit <- rep(NA, length(doseSeq))
      }
      res <- list(fit = mn, se.fit =  as.vector(seFit),
                  residual.scale = sig, df = df)
    }
    res
  }

# Calculate information criteria for nls and lm objects
AIC.fitMCPMod <-
  function(object, ..., k = 2)
  {
    # object - an object of class fitMCPMod
    # k - penalty term
    fm <- object$fm
    if (is.null(fm$weights)) {
      ## user defined model or non-averaged modeling
      return(AIC(fm, k = k))
    }
    RSS <- attr(object, "ResidSS")
    n <- sum(fm$weights)
    sig2 <- RSS/n
    logL <- -n/2*(log(2*pi) + 1 + log(sig2))
    -2*logL + k*(length(coef(object)) + 1) # "+ 1" because of sigma parameter
  }

## model selection
modelSelect <-
  function(data, namSigMod, selMethod, pW, resp, dose, start, nlsControl,
           off, scal, uModPars, addArgs)
  {
    fm <- NA
    warn <- NULL
    nSigMod <- length(namSigMod)
    if (selMethod == "maxT") {  # first maxT option
      i <- 1
      while((length(fm) == 1) && (i <= nSigMod)) {
        nam <- namSigMod[i]
        fm <- fitModel(data, nam, resp, dose, start[[nam]], nlsControl,
                       off, scal, uModPars[[nam]], addArgs[[nam]])
        if(length(fm) == 1) # model didn't converge
          warning(nam, " model did not converge\n")
        i <- i + 1
      }
      if (length(fm) > 1) {  # model converged
        fm <- list(fit = list(fm), base = fm$base)
        attr(fm$fit, "model2") <- names(fm$fit) <- nam
      } else {
        fm <- list(fit = NA, base = NA) # no model converged
      }
    } else {                    # AIC, BIC, aveAIC or aveBIC
      fm <- vector("list", nSigMod)
      crit <- fmBase <- rep(NA, nSigMod)
      if (selMethod == "AIC"| selMethod == "aveAIC") {
        pen <- 2
      } else {
        pen <- log(dim(data)[[1]])
      }
      names(fm) <- names(crit) <- names(fmBase) <- namSigMod
      for(i in 1:nSigMod) {
        nam <- namSigMod[i]
        fitmod <- fitModel(data, nam, resp, dose, start[[nam]], nlsControl,
                           off, scal, uModPars[[nam]], addArgs[[nam]])
        if(!is.list(fitmod)) { # model didn't converge
          fm[[i]] <-  NA
          warning(nam, " model did not converge\n")
        } else { # model converged
          fm[[i]] <- fitmod
          fmBase[i] <- fitmod$base
          crit[i] <- AIC(fitmod, k = pen)
        }
      }
      if (all(is.na(crit))) {
        fm <- NA
      } else {
        if (selMethod == "AIC" | selMethod == "BIC") {
          model2 <- namSigMod[which(crit == min(crit, na.rm = TRUE))]
          fm <- list(fm[[model2]])
          fmBase <- fmBase[model2]
          attr(fm, "model2") <- names(fm) <- model2
          attr(fm, "IC") <- crit
        }
        else { # model averaging
          attr(fm, "model2") <- namSigMod[!is.na(fmBase)]
          attr(fm, "IC") <- crit
          crit <- crit[!is.na(crit)]
          ## subtract const from crit values to avoid numerically 0
          ## values (exp(-0.5*1500)=0!)
          const <- mean(crit)
          if(is.null(pW)){
            pW <- rep(1, length(crit)) # standard 'noninformative' prior
            names(pW) <- names(crit)
          } else {
            pW <- pW[names(crit)]
            if(any(is.na(pW))) stop("pW needs to be a named vector with names equal to the models in the candidate set")
          }
          attr(fm, "weights") <-
            pW*exp(-0.5*(crit-const))/sum(pW*exp(-0.5*(crit-const)))
          attr(fm, "pweights") <- pW
        }
      }
      fm <- list(fit = fm, base = fmBase)
    }
    fm
  }

## dose estimation
getDose <-
  function(dose, ind)
  {
    aa <- !is.na(ind)
    if (!all(aa)) {
      ind <- ind[aa]; dose <- dose[aa]
    }
    if (any(ind)) min(dose[ind])
    else NA
  }

getDoseEst1 <-
  function(fmb,
           CT1,
           false.go.CT1, FGR.CT1,
           false.nogo.CT1, FNGR.CT1,
           CT2,
           false.go.CT2, FGR.CT2,
           false.nogo.CT2 , FNGR.CT2,
           rangeDose, doseEst, selMethod,
           lenDose = 101, uGrad = NULL)
  {
    ## equally spaced points within dose range for prediction
    doseSeq <- seq(rangeDose[1], rangeDose[2], len = lenDose)
    off <- attr(fmb, "off"); scal <- attr(fmb, "scal")
    newdata <- list(dose = doseSeq, off = rep(off, lenDose), scal = rep(scal, lenDose))
    ldePar=1
    val <- double(ldePar)
    base <- fmb$base
    tDose <- matrix(ncol = length(base)-sum(is.na(base)), nrow = 1)
    z <- 1
    for(m in which(!is.na(base))){ # only predict models that converged
      pred <- predict(fmb$fit[[m]], doseSeq, uGrad = uGrad)
      fo   <- data.frame(pred = pred$fit - base[m], se.fit = pred$se.fit)
      df   <- pred$df
      ## selects the desired dose according to MED method

      if(is.na(CT1)){CT1=-Inf}
      if(is.na(CT2)){CT2=-Inf}
      if(false.go.CT1==FALSE){LLCT1=Inf}else{
        gocrt1  <- qt(1-FGR.CT1, df)
        LLCT1<-fo$pred - gocrt1 * fo$se.fit
      }

      if(false.go.CT2==FALSE){LLCT2=Inf}else{
        gocrt2  <- qt(1-FGR.CT2, df)
        LLCT2<-fo$pred - gocrt2 * fo$se.fit
      }

      if(false.nogo.CT1==FALSE){ULCT1=Inf}else{
        nogocrt1  <- qt(1-FNGR.CT1, df)
        ULCT1<-fo$pred + nogocrt1 * fo$se.fit
      }

      if(false.nogo.CT2==FALSE){ULCT2=Inf}else{
        nogocrt2  <- qt(1-FNGR.CT2, df)
        ULCT2<-fo$pred + nogocrt2 * fo$se.fit
      }

      switch(doseEst,
             #"MED1" = ind <- LL >= 0 & UL >= clinRel,
             "MED2" = ind <- LLCT1 >= CT1 & LLCT2 >= CT2,
             #"MED2" = ind <- LL >= 0.15 & fo$pred >= clinRel,
             #"MED3" = ind <- LL <= threshold | fo$pred <= clinRel##mP
             #"MED3" = ind <- LL >= clinRel
             "MED3" = ind <- ULCT1 >= CT1 | ULCT2 >= CT2
             #"MED4" = ind <- UL <= clinRel#MP
      )
      val <- getDose(doseSeq, ind)


      tDose[,z] <- val
      z <- z + 1
    }
    if(is.element(selMethod, c("aveAIC", "aveBIC"))){
      doseAve <- as.vector(tDose%*%attr(fmb$fit, "weights"))
      attr(doseAve, "tdModels") <- tDose
      tDose <- doseAve
    } else {
      tDose <- as.vector(tDose)
    }
    tDose
  }

recovNames <-
  function(names)
  {
    ## function to recover model names (in case of multiple models from one class)
    ## just for use in MCPMod function
    ## example: recovNames(c("emax1", "betaMod", "emax2", "logistic", "usermodel"))
    ## returns: c("emax","betaMod","logistic","usermodel")

    builtIn <- c("linlog", "linear", "quadratic", "emax",
                 "exponential","logistic", "betaMod", "sigEmax")
    newnames <- character()
    i <- 1
    for(nam in names){
      pm <- pmatch(builtIn, nam)
      if(any(!is.na(pm))) {
        newnames[i] <- builtIn[which(!is.na(pm))]
        i <- i+1
      }
      if(all(is.na(pm))) {
        newnames[i] <- nam
        i <- i+1
      }
    }
    unique(newnames)
  }

## function to get reasonable defaults for, off, scal and dePar
getDef <-
  function(off = NULL, scal = NULL, doses, doseEst)
  {
    mD <- max(doses)
    if(is.null(scal)){ # default for scal parameter
      scal <- 1.2*mD
    } else { # check if valid scal is provided
      if(scal < mD){
        stop("'scal' should be >= maximum dose")
      }
    }
    if(is.null(off)){ # default for off parameter
      off <- 0.1*mD
    }
    list(scal = scal, off = off)
  }


### main function implementing methodology, combining several others
MCPMod1 <-
  function(data, models = NULL, contMat = NULL, critV = NULL,
           resp = "resp", dose = "dose", off = NULL, scal = NULL,
           alpha = 0.025, twoSide = FALSE,
           selModel = c("maxT", "AIC", "BIC", "aveAIC", "aveBIC"),
           doseEst = c("MED2", "MED3"),
           std = TRUE, start = NULL,
           uModPars = NULL, addArgs = NULL,
           CT1,
           false.go.CT1, FGR.CT1,
           false.nogo.CT1, FNGR.CT1,
           CT2,
           false.go.CT2, FGR.CT2,
           false.nogo.CT2 , FNGR.CT2, lenDose = 101, pW = NULL,
           control = list(maxiter = 100, tol = 1e-6, minFactor = 1/1024),
           signTtest = 1, pVal = FALSE, testOnly = FALSE,
           mvtcontrol = mvtnorm.control(), na.action = na.fail,
           uGrad = NULL)
  {
    ## data - data.frame with, at least, columns "resp" and "dose"
    ## models - list with component names specifying models
    ##          and optionally elements being used to calculate
    ##          model contrasts; need to have named elements with
    ##          corresponding models and be consistent with contMat,
    ##          if that is non-null
    ## contMat - contrast matrix for candidate set and doses
    ## critV - critical value for MCP test
    ## alpha - significance level for calculating critVal (if = NULL)
    ## selModel - model selection criterion
    ## doseEst - type of dose estimator to be used

    ## pW - named vector of prior probs. for different models

    if (any(is.na(match(c(resp, dose), names(data))))) {
      stop(resp," and/or ", dose, " not found in data")
    }
    data <- na.action(data)
    ind <- match(dose, names(data))
    data <- data[order(data[, ind]), ]
    ## target dose estimate type
    doseEst <- match.arg(doseEst)
    n <- as.vector(table(data[, dose]))           # sample sizes per group
    doses <- sort(unique(data[, dose]))
    ## getting defaults which depend on the DRdata
    def <- getDef(off, scal, doses, doseEst)
    scal <- def[[1]]; off <- def[[2]];
    if (is.null(contMat)) {
      ## need to calculate it
      mu <- modelMeans(models, doses, std, off, scal)
      contMat <- modContr(mu, n)
    }
    ## MCP test first
    tStat <- signTtest * getTstat(data, n, contMat, resp, dose)
    if(twoSide){
      tStat <- abs(tStat)
    }
    if (is.null(critV)) {
      if(pVal) {
        pVals <- pValues(contMat, n, alpha, tStat, mvtcontrol, twoSide)
      }
      critV <- critVal(contMat, n, alpha, mvtcontrol, twoSide = twoSide)
      attr(critV, "Calc") <- TRUE # determines whether cVal was calculated
    } else {
      pVal <- FALSE # pvals are not calculated if critV is supplied
      attr(critV, "Calc") <- FALSE
    }
    indStat <- tStat > critV
    if (!any(indStat) | testOnly) {
      ## only mcp-test or no significant t-stats
      result <- list(signf = any(indStat), model1 = NA, model2 = NA)
    } else {
      ## model selection method
      selMethod <- match.arg(selModel)
      control <- do.call("nls.control", control)
      ## at least one significant, select a model if possible
      namMod <- names(tStat)
      maxTstat <- max(tStat)
      model1 <- namMod[which(tStat == maxTstat)] # model with most sig contrast
      ## significant models, in descending order of tstat
      indSigMod <- 1:length(tStat)
      indSigMod <- indSigMod[rev(order(tStat))][1:sum(indStat)]
      namSigMod <- namMod[indSigMod]       # significant models
      namSigMod <- recovNames(namSigMod)   # remove model nrs.
      # (in case of multiple models)
      fmb <-
        modelSelect(data, namSigMod, selMethod, pW, resp, dose, start,
                    control, off, scal, uModPars, addArgs)
      if (all(is.na(fmb$base))) { # none of sign. model converged
        result <- list(signf = TRUE, model1 = model1, model2 = NA)
      } else {
        ## dose estimation model(s) obtained
        ## move to final step, dose estimation
        ## dose range
        rgDose <- range(doses)
        tDose <- getDoseEst1(fmb, CT1,
                             false.go.CT1, FGR.CT1,
                             false.nogo.CT1, FNGR.CT1,
                             CT2,
                             false.go.CT2, FGR.CT2,
                             false.nogo.CT2 , FNGR.CT2,
                             rgDose, doseEst, selMethod,
                             lenDose, uGrad)
        result <- list(signf = TRUE, model1 = model1,
                       model2 = attr(fmb$fit, "model2"))
      }
    }
    ## add information to the object
    result$input <- list(models=models, resp=resp, dose=dose, off=off,
                         scal=scal, alpha=alpha, twoSide=twoSide,
                         selModel=selModel[1], doseEst=doseEst,
                         std=std, CT1=CT1,
                         false.go.CT1=false.go.CT1,
                         FGR.CT1= FGR.CT1,
                         false.nogo.CT1=false.nogo.CT1,
                         FNGR.CT1=FNGR.CT1,
                         CT2=CT2,
                         false.go.CT2=false.go.CT2,
                         FGR.CT2=FGR.CT2,
                         false.nogo.CT2=false.nogo.CT2 ,
                         FNGR.CT2=FNGR.CT2, uModPars=uModPars,
                         addArgs=addArgs, start = start, uGrad=uGrad,
                         lenDose=lenDose, signTtest=signTtest,
                         pVal=pVal, testOnly=testOnly)
    result$data <- data
    result$contMat <- contMat
    result$corMat <- t(contMat)%*%(contMat/n) /
      sqrt(crossprod(t(colSums(contMat^2/n))))
    result$cVal <- critV
    result$tStat <- tStat
    if(pVal) attr(result$tStat, "pVal") <- pVals
    if (!all(is.na(result$model2))){
      result$fm <- fmb$fit
      result$tdose <- tDose
    }
    oldClass(result) <- "MCPMod"
    result
  }

# print method for MCPMod objects





#### example code
###
### mods <- list(linear = NULL, linlog = NULL, emax = 0.3,
###          quadratic = -0.001, logistic = c(0.4, 0.09))
###
### mn <- c(0.1, 0.4, 0.55, 0.75, 0.9, 1)
### ds <- c(0, 0.05, 0.2, 0.4, 0.6, 1)
### DRdata <- genDFdata(mu = mn, sigma = 0.8, n = 20, doses = ds)
### MM <- MCPMod(DRdata, mods, clinRel = 0.5*max(mn), selModel="aveAIC")
### MM
### summary(MM)
### plot(MM)
###
### user defined model
### emx2 <- function(dose, a, b, d){
###   sigEmax(dose, a, b, d, 1)
### }
### emx2Grad <- function(dose, a, b, d) cbind(1, dose/(dose+d), -b*dose/(dose+d)^2)
### models <- list(emx2=c(0,1,0.1), linear = NULL)
### dats <- genDFdata(mu = mn, doses = ds, sigma=0.8, n=50)
### start <- list(emx2=c(a=0.2, b=0.6, d=0.2))
### MM1 <- MCPMod(dats, models, clinRel = 1, scal = 1.2, selModel="AIC", start = start,
###          uGrad = emx2Grad)

## function to generate DF data
genDFdata <-
  function(model, argsMod, doses, n, sigma, mu = NULL)
  {
    ## generates data.frame with resp and dose columns
    ## corresponding to specified design, mean (either
    ## determined by model or mu) and sigma
    ##
    ## model - string with model function name (first argument
    ##         must be dose, must be vectorized)
    ## argsMod - a named list with the arguments for the model function
    ##
    ## doses - set of doses to be used in the trial
    ## n - sample size (scalar is repeated for each dose)
    ## sigma - error std. deviation
    ## mu - vector of mean values can alternatively specified

    nD <- length(doses)
    dose <- sort(doses)
    if (length(n) == 1) n <- rep(n, nD)
    dose <- rep(dose,  n)
    if(!missing(model)){
      args <- c(list(dose), argsMod)
      mu <- do.call(model, args)
    } else if(!is.null(mu)){
      if(length(doses) != length(mu)){
        stop("'mu' needs to be of the same length as doses")
      }
      mu <- rep(mu,  n)
    } else {
      stop("either 'model' or 'mu' needs to be specified")
    }
    data.frame(dose = dose,
               resp = mu + rnorm(sum(n), sd = sigma))
  }
### examples
### genDFdata("emax", c(e0 = 0.2, eMax = 1, ed50 = 0.05), c(0,0.05,0.2,0.6,1), 20, 1)
### genDFdata(mu = 1:5, doses = 0:4, n = c(20, 20, 10, 5, 1), sigma = 1)


###revised function####
getDoseEst2 <-
  function(fmb,
           CT1.go,
           false.go.CT1, FGR.CT1,
           CT1.nogo,
           false.nogo.CT1, FNGR.CT1,
           CT2.go,
           false.go.CT2, FGR.CT2,
           CT2.nogo,
           false.nogo.CT2 , FNGR.CT2,
           logic.go,
           logic.nogo,
           rangeDose, doseEst, selMethod,
           lenDose = 101, uGrad = NULL,direction='Greater')
  {
    ## equally spaced points within dose range for prediction
    doseSeq <- seq(rangeDose[1], rangeDose[2], len = lenDose)
    off <- attr(fmb, "off"); scal <- attr(fmb, "scal")
    newdata <- list(dose = doseSeq, off = rep(off, lenDose), scal = rep(scal, lenDose))
    ldePar=1
    valgo <- double(ldePar)
    valnogo <- double(ldePar)
    base <- fmb$base
    tDose <- matrix(ncol = length(base)-sum(is.na(base)), nrow = 2)
    z <- 1
    for(m in which(!is.na(base))){ # only predict models that converged
      if(doseEst[1] != "ED"){ # MED estimate
        pred <- predict(fmb$fit[[m]], doseSeq, uGrad = uGrad)
        fo   <- data.frame(pred = pred$fit - base[m], se.fit = pred$se.fit)
        df   <- pred$df
        ## selects the desired dose according to MED method

        if(is.na(CT1.go)){false.go.CT1=FALSE}
        if(is.na(CT2.go)){false.go.CT2=FALSE}
        if(is.na(CT1.nogo)){false.nogo.CT1=FALSE}
        if(is.na(CT2.nogo)){false.nogo.CT2=FALSE}
        if(direction=='Greater'){
          if(false.go.CT1==TRUE){
            gocrt1  <- qt(1-FGR.CT1, df)
            LLCT1<-fo$pred - gocrt1 * fo$se.fit
          }

          if(false.go.CT2==TRUE){
            gocrt2  <- qt(1-FGR.CT2, df)
            LLCT2<-fo$pred - gocrt2 * fo$se.fit
          }

          if(false.nogo.CT1==TRUE){
            nogocrt1  <- qt(1-FNGR.CT1, df)
            ULCT1<-fo$pred + nogocrt1 * fo$se.fit
          }

          if(false.nogo.CT2==TRUE){
            nogocrt2  <- qt(1-FNGR.CT2, df)
            ULCT2<-fo$pred + nogocrt2 * fo$se.fit
          }

          if(doseEst=='APP'){

            if(false.go.CT1==TRUE&false.go.CT2==TRUE){
              if(logic.go=='and'){
                ind1<- LLCT1 >= CT1.go & LLCT2 >= CT2.go
              }
              if(logic.go=='or'){
                ind1<- LLCT1 >= CT1.go | LLCT2 >= CT2.go
              }

            }
            if(false.go.CT1==TRUE&false.go.CT2==FALSE){
              ind1<- LLCT1 >= CT1.go
            }
            if(false.go.CT1==FALSE&false.go.CT2==TRUE){
              ind1<- LLCT2 >= CT2.go
            }

            if(false.nogo.CT1==TRUE&false.nogo.CT2==TRUE){
              if(logic.nogo=='and'){
                ind2<- ULCT1 > CT1.nogo | ULCT2 > CT2.nogo
              }
              if(logic.nogo=='or'){
                ind2<- ULCT1 > CT1.nogo & ULCT2 > CT2.nogo
              }
            }
            if(false.nogo.CT1==TRUE&false.nogo.CT2==FALSE){
              ind2<- ULCT1 > CT1.nogo
            }
            if(false.nogo.CT1==FALSE&false.nogo.CT2==TRUE){
              ind2<- ULCT2 > CT2.nogo
            }










          }
        }
        if(direction=='Less'){
          if(false.go.CT1==TRUE){
            gocrt1  <- qt(1-FGR.CT1, df)
            ULCT1<-fo$pred + gocrt1 * fo$se.fit
          }

          if(false.go.CT2==TRUE){
            gocrt2  <- qt(1-FGR.CT2, df)
            ULCT2<-fo$pred + gocrt2 * fo$se.fit
          }

          if(false.nogo.CT1==TRUE){
            nogocrt1  <- qt(1-FNGR.CT1, df)
            LLCT1<-fo$pred - nogocrt1 * fo$se.fit
          }

          if(false.nogo.CT2==TRUE){
            nogocrt2  <- qt(1-FNGR.CT2, df)
            LLCT2<-fo$pred - nogocrt2 * fo$se.fit
          }

          if(doseEst=='APP'){ # need revision: replace LLCT1<=CT1.go to ULCT1<=CT1.go

            if(false.go.CT1==TRUE&false.go.CT2==TRUE){
              if(logic.go=='and'){
                ind1<- LLCT1 <= CT1.go & LLCT2 <= CT2.go
              }
              if(logic.go=='or'){
                ind1<- LLCT1 <= CT1.go | LLCT2 <= CT2.go
              }

            }
            if(false.go.CT1==TRUE&false.go.CT2==FALSE){
              ind1<- LLCT1 <= CT1.go
            }
            if(false.go.CT1==FALSE&false.go.CT2==TRUE){
              ind1<- LLCT2 <= CT2.go
            }

            if(false.nogo.CT1==TRUE&false.nogo.CT2==TRUE){
              if(logic.nogo=='and'){
                ind2<- ULCT1 < CT1.nogo | ULCT2 < CT2.nogo
              }
              if(logic.nogo=='or'){
                ind2<- ULCT1 < CT1.nogo & ULCT2 < CT2.nogo
              }
            }
            if(false.nogo.CT1==TRUE&false.nogo.CT2==FALSE){
              ind2<- ULCT1 < CT1.nogo
            }
            if(false.nogo.CT1==FALSE&false.nogo.CT2==TRUE){
              ind2<- ULCT2 < CT2.nogo
            }


          }

        }

        valgo <- getDose(doseSeq, ind1)
        valnogo<-getDose(doseSeq,ind2)
      }
      tDose[1,z] <- valgo
      tDose[2,z] <- valnogo
      z <- z + 1
    }
    if(is.element(selMethod, c("aveAIC", "aveBIC"))){
      doseAve <- as.vector(tDose%*%attr(fmb$fit, "weights"))
      attr(doseAve, "tdModels") <- tDose
      tDose <- doseAve
    } else {
      tDose <- tDose
    }
    tDose
  }


MCPMod2 <-
  function(data, models = NULL, contMat = NULL, critV = NULL,
           resp = "resp", dose = "dose", off = NULL, scal = NULL,
           alpha = 0.025, twoSide = FALSE,
           selModel = c("maxT", "AIC", "BIC", "aveAIC", "aveBIC"),
           doseEst = c("MED2", "MED3",'APP'),
           std = TRUE, start = NULL,
           uModPars = NULL, addArgs = NULL,
           CT1.go,
           false.go.CT1, FGR.CT1,
           CT1.nogo,
           false.nogo.CT1, FNGR.CT1,
           CT2.go,
           false.go.CT2, FGR.CT2,
           CT2.nogo,
           false.nogo.CT2 , FNGR.CT2,
           logic.go,logic.nogo,
           lenDose = 101, pW = NULL,
           control = list(maxiter = 100, tol = 1e-6, minFactor = 1/1024),
           signTtest = 1, pVal = FALSE, testOnly = FALSE,
           mvtcontrol = mvtnorm.control(), na.action = na.fail,
           uGrad = NULL,direction='Greater')
  {
    ## data - data.frame with, at least, columns "resp" and "dose"
    ## models - list with component names specifying models
    ##          and optionally elements being used to calculate
    ##          model contrasts; need to have named elements with
    ##          corresponding models and be consistent with contMat,
    ##          if that is non-null
    ## contMat - contrast matrix for candidate set and doses
    ## critV - critical value for MCP test
    ## alpha - significance level for calculating critVal (if = NULL)
    ## selModel - model selection criterion
    ## doseEst - type of dose estimator to be used

    ## pW - named vector of prior probs. for different models

    if (any(is.na(match(c(resp, dose), names(data))))) {
      stop(resp," and/or ", dose, " not found in data")
    }
    data <- na.action(data)
    ind <- match(dose, names(data))
    data <- data[order(data[, ind]), ]
    ## target dose estimate type
    doseEst <- match.arg(doseEst)
    n <- as.vector(table(data[, dose]))           # sample sizes per group
    doses <- sort(unique(data[, dose]))
    ## getting defaults which depend on the DRdata
    def <- getDef(off, scal, doses, doseEst)
    scal <- def[[1]]; off <- def[[2]];
    if (is.null(contMat)) {
      ## need to calculate it
      mu <- modelMeans(models, doses, std, off, scal)
      contMat <- modContr(mu, n)
    }
    ## MCP test first
    tStat <- signTtest * getTstat(data, n, contMat, resp, dose)
    if(twoSide){
      tStat <- abs(tStat)
    }
    if (is.null(critV)) {
      if(pVal) {
        pVals <- pValues(contMat, n, alpha, tStat, mvtcontrol, twoSide)
      }
      critV <- critVal(contMat, n, alpha, mvtcontrol, twoSide = twoSide)
      attr(critV, "Calc") <- TRUE # determines whether cVal was calculated
    } else {
      pVal <- FALSE # pvals are not calculated if critV is supplied
      attr(critV, "Calc") <- FALSE
    }
    indStat <- tStat > critV
    if (!any(indStat) | testOnly) {
      ## only mcp-test or no significant t-stats
      result <- list(signf = any(indStat), model1 = NA, model2 = NA)
    } else {
      ## model selection method
      selMethod <- match.arg(selModel)
      control <- do.call("nls.control", control)
      ## at least one significant, select a model if possible
      namMod <- names(tStat)
      maxTstat <- max(tStat)
      model1 <- namMod[which(tStat == maxTstat)] # model with most sig contrast
      ## significant models, in descending order of tstat
      indSigMod <- 1:length(tStat)
      indSigMod <- indSigMod[rev(order(tStat))][1:sum(indStat)]
      namSigMod <- namMod[indSigMod]       # significant models
      namSigMod <- recovNames(namSigMod)   # remove model nrs.
      # (in case of multiple models)
      fmb <-
        modelSelect(data, namSigMod, selMethod, pW, resp, dose, start,
                    control, off, scal, uModPars, addArgs)
      if (all(is.na(fmb$base))) { # none of sign. model converged
        result <- list(signf = TRUE, model1 = model1, model2 = NA)
      } else {
        ## dose estimation model(s) obtained
        ## move to final step, dose estimation
        ## dose range
        rgDose <- range(doses)
        tDose <- getDoseEst2(fmb, CT1.go,
                             false.go.CT1, FGR.CT1,
                             CT1.nogo,
                             false.nogo.CT1, FNGR.CT1,
                             CT2.go,
                             false.go.CT2, FGR.CT2,
                             CT2.nogo,
                             false.nogo.CT2 , FNGR.CT2,
                             logic.go,logic.nogo,
                             rgDose, doseEst, selMethod,
                             lenDose, uGrad,direction)
        result <- list(signf = TRUE, model1 = model1,
                       model2 = attr(fmb$fit, "model2"))
      }
    }
    ## add information to the object
    result$input <- list(models=models, resp=resp, dose=dose, off=off,
                         scal=scal, alpha=alpha, twoSide=twoSide,
                         selModel=selModel[1], doseEst=doseEst,
                         std=std, CT1.go=CT1.go,
                         false.go.CT1=false.go.CT1,
                         FGR.CT1= FGR.CT1,
                         CT1.nogo=CT1.nogo,
                         false.nogo.CT1=false.nogo.CT1,
                         FNGR.CT1=FNGR.CT1,
                         CT2.go=CT2.go,
                         false.go.CT2=false.go.CT2,
                         FGR.CT2=FGR.CT2,
                         CT2.nogo=CT2.nogo,
                         false.nogo.CT2=false.nogo.CT2 ,
                         FNGR.CT2=FNGR.CT2,
                         logic.go=logic.go,logic.nogo=logic.nogo,
                         uModPars=uModPars,
                         addArgs=addArgs, start = start, uGrad=uGrad,
                         lenDose=lenDose, signTtest=signTtest,
                         pVal=pVal, testOnly=testOnly)
    result$data <- data
    result$contMat <- contMat
    result$corMat <- t(contMat)%*%(contMat/n) /
      sqrt(crossprod(t(colSums(contMat^2/n))))
    result$cVal <- critV
    result$tStat <- tStat
    if(pVal) attr(result$tStat, "pVal") <- pVals
    if (!all(is.na(result$model2))){
      result$fm <- fmb$fit
      result$tdose <- tDose
    }
    oldClass(result) <- "MCPMod"
    result
  }


MCPMod_MED_APP_parallel = function(datatype="Binomial",dose, case,
                                   mean,sd,
                                   mods,
                                   n,
                                   alpha=0.1,
                                   CT1.go=0,
                                   false.go.CT1=TRUE, FGR.CT1=0.1,
                                   CT1.nogo=0.1,
                                   false.nogo.CT1=TRUE, FNGR.CT1=0.3,
                                   CT2.go=0.15,
                                   false.go.CT2=TRUE, FGR.CT2=0.5,
                                   CT2.nogo=0.15,
                                   false.nogo.CT2=TRUE, FNGR.CT2=0.9,
                                   logic.go='and',logic.nogo='or',
                                   overlap.option='GO',
                                   nsim=10,seednum=369,direction='Greater'){


  ######### function #####################################################


  len<-length(dose)
  if(len!=length(n)){
    stop('The number of dose levels and the length of sample size in each dose level are different')
  }
  datadose<-rep(dose,n)

  ncore <- detectCores()
  cl<-makeCluster(ncore)
  registerDoParallel(cl)
  results<-foreach(i = 1:nsim,.packages=c('MCPAPP','mvtnorm',
                                          'lattice'),.combine=rbind) %dopar% {
                                            index=0
                                            set.seed(seednum*i)
                                            dataresp<-rep(NA,length(datadose))
                                            if(datatype=="Binomial"){
                                              for(ii in 1:len){
                                                dataresp[(index+1):(index+n[ii])]<-rbinom(n[ii],1,case[ii])
                                                index=index+n[ii]
                                              }
                                            }
                                            if(datatype=='Normal'){
                                              for(ii in 1:len){
                                                dataresp[(index+1):(index+n[ii])]<-rnorm(n[ii],mean=mean[ii],sd=sd[ii])
                                                index=index+n[ii]
                                              }
                                            }
                                            data = data.frame(dose=datadose, resp=dataresp)

                                            df = MCPMod2(data, mods, alpha=alpha,
                                                         CT1.go=CT1.go,
                                                         false.go.CT1=false.go.CT1, FGR.CT1=FGR.CT1,
                                                         CT1.nogo=CT1.nogo,
                                                         false.nogo.CT1=false.nogo.CT1, FNGR.CT1=FNGR.CT1,
                                                         CT2.go=CT2.go,
                                                         false.go.CT2=false.go.CT2, FGR.CT2=FGR.CT2,
                                                         CT2.nogo=CT2.nogo,
                                                         false.nogo.CT2=false.nogo.CT2 , FNGR.CT2=FNGR.CT2,
                                                         logic.go=logic.go,logic.nogo=logic.nogo,
                                                         pVal=FALSE,
                                                         selModel="maxT",doseEst="APP",  scal = NULL,direction=direction)

                                            if(is.null(df$tdose)){
                                              minED = -999
                                              minED2=-999

                                            }else{
                                              minED=df$tdose[1]
                                              minED2=df$tdose[2]}
                                            c(minED,minED2)
                                          }
  stopCluster(cl)
  target_dose <-results[,1]
  target_dose2<-results[,2]

  target_dose[is.na(target_dose)]=-1000
  target_dose2[is.na(target_dose2)]=-1000

  go_index=1-(target_dose==-1000|target_dose==-999)
  nogo_index=1*(target_dose2==-1000|target_dose2==-999)
  notrend_index=1*(target_dose==-999)
  t=go_index+nogo_index
  overlap.flag=1*(sum(t==2)>0)
  if(overlap.option=='GO'){
    if(sum(go_index==1)>0){
      nogo_index[go_index==1]<-0
    }
  }
  if(overlap.option=='NOGO'){
    if(sum(nogo_index==1)>0){
      go_index[nogo_index==1]<-0
    }
  }

  go=c(sum(go_index, na.rm = TRUE))/nsim
  nogo=c(sum(nogo_index, na.rm = TRUE))/nsim
  grey=1-go-nogo
  notrend=c(sum(notrend_index,na.rm=TRUE))/nsim
  output = list(n=n, go=go, nogo=nogo, grey=grey,notrend=notrend,
                overlap.flag=overlap.flag,overlap.option=overlap.option)
  return(output)
}




MCPMod_MED_APP_interim_parallel = function(kkk,datatype="Binomial",dose, case,
                                           mean,sd,
                                           mods,
                                           n,
                                           num_interim=3,
                                           alpha=0.1,
                                           CT1.go=c(0,0,0),
                                           false.go.CT1=c(TRUE,TRUE,TRUE), FGR.CT1=c(0.1,0.1,0.1),
                                           CT1.nogo=c(0.1,0.1,0.1),
                                           false.nogo.CT1=c(TRUE,TRUE,TRUE), FNGR.CT1=c(0.3,0.3,0.3),
                                           CT2.go=c(0.15,0.15,0.15),
                                           false.go.CT2=c(TRUE,TRUE,TRUE), FGR.CT2=c(0.5,0.5,0.5),
                                           CT2.nogo=c(0.15,0.15,0.15),
                                           false.nogo.CT2=c(TRUE,TRUE,TRUE), FNGR.CT2=c(0.9,0.9,0.9),
                                           logic.go=c('and','and','or'),logic.nogo=c('or','and','and'),
                                           task=c('Futility','Superiority','Futility and superiority'),
                                           overlap.option=c('GO','NOGO','GO'),
                                           nsim_IA=10,seednum=369,
                                           direction=c('Greater','Less','Greater')){


  ######### function #####################################################


  len<-length(dose)
  if(len!=length(n[num_interim,])){
    stop('The number of dose levels and the length of sample size in each dose level are different')
  }


  go_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
  nogo_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
  inconclusive_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
  t<-matrix(NA,nrow=nsim_IA,ncol=num_interim)
  IA_go_matrix<-matrix(NA,nrow=nsim_IA,ncol=num_interim) ###whether continue to next stage
  overlap=rep(NA,num_interim)
  temptable=c()
  table<-matrix(NA,ncol=num_interim+1,nrow=7)
  for(interim_index in 1:num_interim){
    ncore <- detectCores()
    cl<-makeCluster(ncore)
    registerDoParallel(cl)
    results<-foreach(i = 1:nsim_IA,.packages=c('MCPAPP','mvtnorm',
                                               'lattice'),.combine=rbind) %dopar% {
                                                 index=0

                                                 datadose<-rep(dose,n[interim_index,])
                                                 dataresp<-rep(NA,length(datadose))

                                                 if(datatype=="Binomial"){
                                                   for(ii in 1:len){
                                                     set.seed(seednum*i*ii)
                                                     dataresp[(index+1):(index+n[interim_index,ii])]<-rbinom(n[interim_index,ii],1,case[ii])
                                                     index=index+n[interim_index,ii]
                                                   }
                                                 }
                                                 if(datatype=='Normal'){
                                                   for(ii in 1:len){
                                                     set.seed(seednum*i*ii)
                                                     dataresp[(index+1):(index+n[interim_index,ii])]<-rnorm(n[interim_index,ii],mean=mean[ii],sd=sd[ii])
                                                     index=index+n[interim_index,ii]
                                                   }
                                                 }

                                                 data = data.frame(dose=datadose, resp=dataresp)

                                                 df = MCPMod2(data, mods, alpha=alpha,
                                                              CT1.go=CT1.go[interim_index],
                                                              false.go.CT1=false.go.CT1[interim_index], FGR.CT1=FGR.CT1[interim_index],
                                                              CT1.nogo=CT1.nogo[interim_index],
                                                              false.nogo.CT1=false.nogo.CT1[interim_index], FNGR.CT1=FNGR.CT1[interim_index],
                                                              CT2.go=CT2.go[interim_index],
                                                              false.go.CT2=false.go.CT2[interim_index], FGR.CT2=FGR.CT2[interim_index],
                                                              CT2.nogo=CT2.nogo[interim_index],
                                                              false.nogo.CT2=false.nogo.CT2[interim_index] , FNGR.CT2=FNGR.CT2[interim_index],
                                                              logic.go=logic.go[interim_index],logic.nogo=logic.nogo[interim_index],
                                                              pVal=FALSE,
                                                              selModel="maxT",doseEst="APP",  scal = NULL,direction=direction[interim_index])

                                                 if(is.null(df$tdose)){
                                                   minED = -999
                                                   minED2=-999

                                                 }else{
                                                   minED=df$tdose[1]
                                                   minED2=df$tdose[2]}
                                                 c(minED,minED2)


                                               }
    stopCluster(cl)
    target_dose <-results[,1]
    target_dose2<-results[,2]
    #print(cbind(target_dose,target_dose2))
    target_dose[is.na(target_dose)]=-1000
    target_dose2[is.na(target_dose2)]=-1000

    go_matrix[,interim_index]=1-(target_dose==-1000|target_dose==-999)
    nogo_matrix[,interim_index]=1*(target_dose2==-1000|target_dose2==-999)
    #notrend_index[,interim_index]=1*(target_dose==-999)
    t[,interim_index]=go_matrix[,interim_index]+nogo_matrix[,interim_index]
    overlap[interim_index]=1*(sum(t[,interim_index]==2)>0)
    if(overlap.option[interim_index]=='GO'){
      if(sum(go_matrix[,interim_index]==1)>0){
        nogo_matrix[go_matrix[,interim_index]==1,interim_index]<-0
      }
    }
    if(overlap.option[interim_index]=='NOGO'){
      if(sum(nogo_matrix[,interim_index]==1)>0){
        go_matrix[nogo_matrix[,interim_index]==1,interim_index]<-0
      }
    }

    inconclusive_matrix[,interim_index]<-rep(1,nsim_IA)-go_matrix[,interim_index]-nogo_matrix[,interim_index]

  }

  for(ii in 1:(num_interim)){
    if(task[ii]=='Futility'){
      IA_go_matrix[,ii]=inconclusive_matrix[,ii]+go_matrix[,ii]
    }
    if(task[ii]=='Superiority'){
      IA_go_matrix[,ii]=inconclusive_matrix[,ii]+nogo_matrix[,ii]
    }
    if(task[ii]=='Futility and superiority'){
      IA_go_matrix[,ii]=inconclusive_matrix[,ii]
    }
  }


  cum_IA_go_matrix<-t(apply(IA_go_matrix,1,cumprod))

  for(j in 1:(num_interim)){
    table[1,j]=paste0(n[j,],collapse = ',')
    if(datatype=='Binomial'){
      table[2,j]=paste0('RR in each dose: ',paste0(case,collapse=','))
    }
    if(datatype=='Normal'){
      table[2,j]=paste0('Mean in each dose: ',paste0(mean,collapse=','),';','sd in each dose: ',paste0(sd,collapse=',',':'))
    }
    table[3,j]=task[j]
    if(j==1){
      if(task[j]=='Superiority'|task[j]=='Futility and superiority'){
        table[4,j]=round(sum(go_matrix[,j]==1)/nsim_IA,3)}else{table[4,j]=0}
      if(task[j]=='Futility'|task[j]=='Futility and superiority'){
        table[6,j]=round(sum(nogo_matrix[,j]==1)/nsim_IA,3)
      }else{table[6,j]=0}
    }else{
      if(task[j]=='Superiority'|task[j]=='Futility and superiority'){
        table[4,j]=round(sum(go_matrix[,j]==1&cum_IA_go_matrix[,j-1]==1)/nsim_IA,3)}else{table[4,j]=0}
      if(task[j]=='Futility'|task[j]=='Futility and superiority'){
        table[6,j]=round(sum(nogo_matrix[,j]==1&cum_IA_go_matrix[,j-1]==1)/nsim_IA,3)
      }else{table[6,j]=0}
    }
    table[5,j]=round(sum(cum_IA_go_matrix[,j]==1)/nsim_IA,3)
    table[7,j]<-ifelse(overlap[j]==1,paste0('GO/NOGO zones overlapped, classified by criterion of ',overlap.option[j]),'None')

  }
  table[1,num_interim+1]=''
  if(datatype=='Binomial'){
    table[2,num_interim+1]=paste0('RR in each dose: ',paste0(case,collapse=','))
  }
  if(datatype=='Normal'){
    table[2,num_interim+1]=paste0('Mean in each dose: ',paste0(mean,collapse=','),';','sd in each dose: ',paste0(sd,collapse=',',':'))
  }
  table[3,num_interim+1]=''
  table[4,num_interim+1]=round(sum(as.numeric(table[4,1:num_interim])),3)
  table[5,num_interim+1]=round(as.numeric(table[5,num_interim]),3)
  table[6,num_interim+1]=round(sum(as.numeric(table[6,1:num_interim])),3)
  table[7,num_interim+1]=''

  table<-as.table(table)
  tablecolname<-c(paste0('Interim analysis ',1:(num_interim-1)),'Final analysis',"Summary")
  tablerowname<-c('Sample size','True parameters','Task','Success','To next interim/final or inconclusive',
                  'Stop','Warning')
  table<-cbind(tablerowname,rep(kkk,7),table)
  colnames(table)<-c(" ",'Setting',tablecolname)
  return(table)
}



