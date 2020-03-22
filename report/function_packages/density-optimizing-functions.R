
## Perform numerical optimization to the density functionand save the result to the given path
optimize.density <- function(X, density, optim.out.pathname, gradient = NULL, control = NULL, echo = T) {

  opt.density <- function(x) {
    density(t(as.matrix(x)))
  }
  
  if (is.null(control)) {
    control <- list("maxit" = 10000000, "fnscale" = -1)
  }
  
  optim.out.list <- NULL
  
  if (file.exists(optim.out.pathname)) {
    source(optim.out.pathname)
  }
  
  i.start <- length(optim.out.list) + 1
  for (i in i.start:nrow(X)) {
    optim.out <- optim(t(as.matrix(X[i,])), fn = opt.density, gr = gradient, method = "BFGS", control = control)
    if (is.null(optim.out.list)) optim.out.list <- list(optim.out)
    else optim.out.list <- c(optim.out.list, list(optim.out))
    if (i %% 10 == 0) dump("optim.out.list", optim.out.pathname)
    if (echo == T) cat(" ", i)
  }
  dump("optim.out.list", optim.out.pathname)
}

## Read optimizing results
read.optim <- function(X, optim.out.pathname) {
  source(optim.out.pathname)
  n <- length(optim.out.list)
  m <- ncol(X)
  
  convergence <- rep(0, n)
  optima <- matrix(0, nrow = n, ncol = m)
  for (i in 1:n) {
    convergence[i] <- optim.out.list[[i]]$convergence
    optima[i,] <- optim.out.list[[i]]$par
  }
  print(sum(convergence))
  return(optima)
}

