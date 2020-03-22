## Cluster Analysis and Other Unsupervised Learning Methods
## STAT 593B, Spring 2014
## Werner Stuetzle

##=================================================================
## Functions for kernel density estimation and cross-validation
##================================================================= 

standardize <- function(X) {
  centered <- sweep(X, 2, apply(X, 2, mean))
  sds <- sqrt(apply(centered, 2, var))
  scaled <- sweep(centered, 2, sds, FUN = "/")
  return(scaled)
}

sphere <- function(X) {
  Centered <- sweep(X, 2, apply(X, 2, mean))
  Sigma  <- var(Centered)
  R <- chol(Sigma)
  Sphered <- t(solve(t(R), t(Centered)))
  return(Sphered)
}

##-----------------------------------------------------------------

make.gaussian.kernel.density.estimate <- function(X.train, h) {
  density.estimate <- function(X.eval) {
    phat <- gaussian.kernel.density.estimate(X.train, X.eval, h, cv = F)
    return(phat)
  }
  return(density.estimate)
}

##-----------------------------------------------------------------

gaussian.kernel.density.estimate <- function(X.obs, X.eval, h, cv = F) {
  ## K <- function(tsquared, h) {
  ##     return(exp(-tsquared / (2 * h^2)) / h^p)
  ##   }
  K <- function(tsquared, h) {
    return(exp(-tsquared / (2 * h^2)) / (((2 * pi)^(p/2)) * h^p))
  }
  if (is.vector(X.obs)) {
    X.obs <- matrix(X.obs, ncol = 1)
    X.eval <- matrix(X.eval, ncol = 1)
  }
  n.eval <- nrow(X.eval)
  n.obs <- nrow(X.obs)
  p <- ncol(X.obs)
  dens <- rep(0, n.eval)
  squared.dist <- gsl.interpoint.distance.matrix(X.eval, X.obs) 
  for (i in 1:n.eval) {
    kernel.values <- K(squared.dist[i,], h)
    dens[i] <- sum(kernel.values) / n.obs
    if (cv) {
      dens[i] <- (sum(kernel.values) - kernel.values[i]) / (n.obs - 1)
    }
  }
  return(dens)
}

##-----------------------------------------------------------------

gaussian.kernel.least.squares.cv <- function(X, h) {
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  n <- nrow(X)
  p <- ncol(X)
  nh <- length(h)
  cv <- rep(0, nh)
  K <- function(tsquared, h) {
    return(exp(-tsquared / (2 * h^2)) / h^p)
  }
  Kstar <- function(tsquared, h) {
    return(K(tsquared, h*sqrt(2)))
  }
  D <- gsl.interpoint.distance.matrix(X, X)
  for (i in 1:nh) {
    cv[i] <- sum(Kstar(D, h[i]))/n^2 - 2 * sum(K(D, h[i]))/(n * (n-1)) +
      2 * K(0, h[i])/(n-1)
  }
  return(cv)
}

##-----------------------------------------------------------------

golden.section <- function(fun, xleft = 0.05, xright = 1, ngrid = 10, maxit = 8){
  xgrid <- seq(xleft, xright, length = ngrid)
  fgrid <- rep(0, ngrid)
  for (i in 1:ngrid) fgrid[i] <- fun(xgrid[i])
  iopt <- which.min(fgrid)
  if ((iopt == 1)|(iopt == ngrid)) return(-1)
  xleft <- xgrid[iopt-1]; fleft = fgrid[iopt-1]
  xmiddle <- xgrid[iopt]; fmiddle = fgrid[iopt]
  xright <- xgrid[iopt+1]; fright = fgrid[iopt+1]
  ## We have initial bracket
  for (i in 1:maxit) {
    if ((xmiddle - xleft) > (xright - xmiddle)) {
      xnew <- xmiddle - 0.38 * (xmiddle - xleft)
      fnew <- fun(xnew)
      if (fnew < fmiddle) {
        xright <- xmiddle; fright <- fmiddle; xmiddle <- xnew; fmiddle <- fnew
      }
      else {
        xleft <- xnew; fleft = fnew
      }
    }
    else {
      xnew <- xright - 0.38 * (xright - xmiddle)
      fnew <- fun(xnew)
      if (fnew < fmiddle) {
        xleft <- xmiddle; fleft <- fmiddle; xmiddle <- xnew; fmiddle <- fnew
      }
      else{
        xright <- xnew; fright <- fnew
      }
    }
  }
  return(xmiddle)
}

##-----------------------------------------------------------------

cv.search <- function(X, trial.par = c(0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4,
                           0.6, 0.8, 1.0, 1.5, 2.0, 4, 8)) {
  cv.function <- gaussian.kernel.least.squares.cv
  search.for.h <- T
  cv <- cv.function(X, trial.par)
  finite <- is.finite(cv)
  opt.smopar <- NA
  finite.par <- NA
  finite.cv <- NA
  if (sum(finite) > 3) {
    finite.cv <- cv[finite]
    finite.par <- trial.par[finite]
    imin <- which.min(finite.cv)
    if((imin != 1) & (imin != length(finite.cv))) {
      opt.smopar <- finite.par[imin]
      if (search.for.h) {
        fun <- function(h) {
          return(cv.function(X, h))
        }
        opt.smopar <- golden.section(fun, xleft = finite.par[imin - 1],
                                     xright = finite.par[imin + 1], maxit = 8)
      }
    }
  }
  return(list(opt.smopar = opt.smopar, smopars = finite.par, cv = finite.cv))
}

##-----------------------------------------------------------------

gsl.interpoint.distance.matrix <- function (X, Y) {
  p <- ncol(X)
  nx <- nrow(X)
  ny <- nrow(Y)
  norm2.x <- apply(X^2, 1, sum)
  norm2.y <- apply(Y^2, 1, sum)
  dist <- matrix(norm2.x, nrow = nx, ncol = ny) +
    matrix(norm2.y, nrow = nx, ncol = ny, byrow = T) -
      2 * X %*% t(Y)
  return(dist)
}


##------------------------------------------------------------------
## Make gradient function of gaussian kernel density estimator
## Added by Richard Zhu
make.kde.gradient <- function (X.train, h) {
  gradient.function <- function(X.eval) {
    gradient <- gaussian.gradient(X.train, X.eval, h)
    return(gradient)
  }
  return(gradient.function)
}

##------------------------------------------------------------------
## Find gradient of gaussian kernel density estimator
## Added by Richard Zhu
gaussian.gradient <- function (X.obs, X.eval, h) {
  if (is.vector(X.obs)) {
    X.obs <- matrix(X.obs, ncol = 1)
    X.eval <- matrix(X.eval, ncol = 1)
  }
  if (is.vector(X.eval)) {
    X.eval <- t(matrix(X.eval, ncol = 1))
  }
  n.eval <- nrow(X.eval)
  n.obs <- nrow(X.obs)
  p <- ncol(X.obs)
  squared.dist <- gsl.interpoint.distance.matrix(X.eval, X.obs) 
  K <- function(tsquared, h) {
    return(exp(-tsquared / (2 * h^2)) / (((2 * pi)^(p/2)) * h^p))
  }
  grad <- rep(0, p)
  for (i in 1:n.eval) {
    kernel.values <- K(squared.dist[i,], h)
    for (j in 1:p) {
      dif <- -(X.eval[i, j] - X.obs[1:n.obs, j])
      grad[j] <- sum(kernel.values * dif / h^2) / n.obs
    }
  }
  return(grad)
}


##-----------------------------------------------------------------
## Example for use of code

## dir.sep <- "/"
## data.dir <- "http://www.stat.washington.edu/wxs/Unsupervised-learning-spring-2016/Data"
## tc.dir <- paste(data.dir, "Test-collection", sep = dir.sep)

## source(paste(tc.dir, "nice-univariate-trimodal.R", sep = dir.sep), echo = echo)

## X.sphered <- sphere(X)
## cv.search.out <- cv.search(X.sphered)
## optimal.bandwidth <- cv.search.out$opt.smopar
## optimal.bandwidth

## grid <- matrix(seq(from = min(X.sphered), to = max(X.sphered),
##                    length.out = 500), ncol = 1)
## phat <- gaussian.kernel.density.estimate(X.sphered, grid, optimal.bandwidth)
## plot(grid, phat, type = "l")


## Alternatively

## gkde <- make.gaussian.kernel.density.estimate(X.sphered, optimal.bandwidth)
## phat <- gkde(grid)
## plot(grid, phat, type = "l")




