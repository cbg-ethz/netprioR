#' Conjugate Gradient Solver
#' 
#' Solves linear equation systems iteratively
#' 
#' @author Fabian Schmich
#' @import Matrix
#' @export
#' 
#' @param A Matrix
#' @param b Coefficients
#' @param x0 Starting solution
#' @param threshold Termination threshold
#' @param verbose Show iterative progress
#' @return Solution for equation system
conjugate.gradient <- function(A, b, x0, threshold = 1e-20, verbose = TRUE) {
  r0 <- A %*% x0 - b
  p0 <- -r0
  
  k <- 0
  rk <- r0
  pk <- p0
  xk <- x0
  
  while (as.numeric(t(rk) %*% rk) > threshold) {
    ak <- as.numeric(t(rk) %*% rk / t(pk) %*% A %*% pk)
    xkp1 <- xk + ak * pk
    rkp1 <- rk + ak * A %*% pk
    bkp1 <- as.numeric(t(rkp1) %*% rkp1 / t(rk) %*% rk)
    pkp1 <- -rkp1 + bkp1 * pk
    
    if (verbose)  print(sprintf("r%d: %.10f", k, as.numeric(t(rk) %*% rk)))
    
    k <- k + 1
    rk <- rkp1
    pk <- pkp1
    xk <- xkp1
  }
  return(xkp1)
}

#' Laplacian Matrix
#' 
#' Compute the laplacian matrix
#' 
#' @author Fabian Schmich
#' @export
#' @import Matrix
#' 
#' @param x Input matrix
#' @return Laplacian matrix
laplacian <- function(x) {
  stopifnot(nrow(x) == ncol(x)) 
  r <- Diagonal(x=(x %*% rep(1, nrow(x)))[,1])
  dimnames(r) <- dimnames(x)
  return(r - x)
}

#' Make regular Grid
#' 
#' Creates a grid of n x m points
#' 
#' @author Fabian Schmich
#' @export
#' 
#' @param n Height of grid
#' @param m Width of grid
#' @param graph return graph instead of adjacency matrix
#' @return Grid graph or adjacency matrix
makegrid <- function(n, m, dimnames = NULL) {
  a <- n * m
  E <- Matrix(0, nrow = a, ncol = a)
  E[cbind(1:(a-1), 2:a)] <- E[cbind(2:a, 1:(a-1))] <- 1
  for (i in seq(2 * n, a, by = n)) {
    d <- cbind(i:(i - 2 * n + 1), (i - 2 * n + 1):i)
    E[d] <- 1
  }
  if (!is.null(dimnames)) dimnames(E) <- dimnames
  return(E)
}
