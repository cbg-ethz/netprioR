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
conjugate_gradient <- function(A, b, x0, threshold = 1e-20, verbose = TRUE) {
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
  if (!is.null(dimnames(x))) dimnames(r) <- dimnames(x)
  return(r - x)
}

#' Compute the bandwidth of a matrix
#' 
#' @author Fabian Schmich
#' @export
#' @import Matrix
#' 
#' @param x Inpute matrix
#' @return Bandwidth
bandwidth <- function(x) {
  k1 <- sapply(1:nrow(x), function(r) r - min(which(x[r,] != 0)))
  k2 <- sapply(1:nrow(x), function(r) max(which(x[r,] != 0)) - r)
  k <- pmax(k1, k2)
  return(max(k))
}

#' Cuthill McKee (CM) algorithm
#' 
#' Transform sparse matrix into a band matrix
#' 
#' @author Fabian Schmich
#' @export
#'  
#' @param x Input matrix
#' @return Band matrix
cuthill_mckee <- function(x) {
  degs <- data.frame(Idx=1:ncol(x), NonZero=apply(x, 1, function(x) length(which(x != 0))))
  R <- degs$Idx[which.min(degs$NonZero)]
  i <- 1
  for (i in 1:ncol(x)) {
    Ai <- setdiff(which(x[R[i],] != 0), R)
    if (length(Ai) > 0) {
      Ai <- Ai[order(degs$NonZero[Ai], decreasing=F)]
      R <- append(R, Ai)
    } else {
      R <- append(R, degs$Idx[-R][which.min(degs$NonZero[-R])])
    }
    i <- i + 1
  }
  rR <- rev(R)
  return(x[rR, rR])
}


#' Create a grid graph of n x m points
#' 
#' @author Fabian Schmich
#' @export
#' @import Matrix
#' 
#' @param n Height of grid
#' @param m Width of grid
#' @param graph return graph instead of adjacency matrix
#' 
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

#' Construct precision matrix for GMRF
#' 
#' @author Fabian Schmich
#' @export
#' @import Matrix
#' 
#' @param V Vector of (named) vertices
#' @param E Matrix of similarities / affinities, i.e. edges between vertices
#' @param labels String sequence of labels for vertices
#' @param theta Parameter to tune label bias
#' 
#' @return Precision matrix Q
pmat <- function(genes, interactions, label.priors, theta) {
  # Number of vertices and classes
  N <- length(genes)
  L <- length(label.priors)
  # Assertions
  stopifnot(all(dim(interactions) == N)) # square matrix with dimensons |nodes|
  stopifnot(!is.null(names(label.priors)))
  stopifnot(is.list(label.priors) & all(lapply(label.priors, length) == N)) # label lists must provide affinity for each node
  stopifnot(all(rowSums(do.call("cBind", label.priors)) == 1)) # label prior over all classes must sum to 1
  # Construct precision matrix
  dn <- NULL
  if (is.null(dn)) dn <- c(paste("N", 1:N, sep = ""), names(label.priors))
  C <- (do.call("cBind", label.priors))
  W <- Matrix(rBind(cBind(interactions, theta * C), cBind(theta * t(C), diag(length(label.priors)))), dimnames = list(dn, dn))
  Q <- laplacian(W)
#   Q <- Matrix(0, nrow = N + L, ncol = N + L, dimnames = list(dn, dn))
#   Q[1:N, 1:N] <- laplacian(E) + theta * diag(N)
#   Q[1:N, (N+1):(N+L)] <- -theta * (do.call("cbind", label.prior))
#   Q[(N+1):(N+L), 1:N] <- t(Q[1:N, (N+1):(N+L)])
#   Q[(N+1):nrow(Q), (N+1):ncol(Q)] <- Matrix(0, nrow = L, ncol = L)
  
  return(Q)
}

#' Sampling from a GMRF
#' 
#' From "GMRF Theory and Application" Alogirhtm 2.5 and Algorithm 2.6
#' 
#' @author Fabian Schmich
#' 
#' @export
#' 
#' @param n Number of samples to be drawn
#' @param Q Precision Matrix
#' @param M Indices of nodes to condition on
#' @param mu Mean of GMRF
#' @param A Constraint: Left side
#' @param e Constraint: Right side
#' 
#' @return n samples from GMRF
gmrf <- function(n, Q, M, xM, mu = 0, constr = list(A, e), use.inla = FALSE) {
  if (length(mu) == 1) mu <- rep(mu, ncol(Q))
  mu <- as.numeric(mu)
  if (missing(M)) { # Unconditional sampling
    Lt <- chol(Q) #1
    L <- t(Lt)    
    z <- rBind(mvrnorm(n = n, mu = mu, Sigma = diag(ncol(Q)))) #2
    v <- apply(z, 1, function(x) backsolve(r = Lt, x = x)) #3    
    if (missing(constr)) { # Algorithm 2.5 (unconstrained)
      if (use.inla) {
        inla.qsample(n = n, Q = Q, mu = mu, logdens =  FALSE)
      } else {
        w <- forwardsolve(l = L, x = mu)
        return(v + backsolve(r = Lt, x = w))
      }      
    } else { # Kriging (constrained)
      if (use.inla) {
        inla.qsample(n = n, Q = Q, mu = mu, constr = constr, logdens =  FALSE)
      } else { 
        x <- mu + v #4
        V <- apply(constr$A, 1, function(x) { #5
          v <- forwardsolve(l = L, x = x)
          backsolve(r = Lt, x = v)
        })
        W <- constr$A %*% V #6
        WLt <- chol(W)
        WL <- t(WLt)
        U <- rBind(apply(V, 1, function(x) { #7
          u <- forwardsolve(l = WL, x = x)
          backsolve(r = WLt, x = u)
        }))
        c <- constr$A %*% x - constr$e #8
        return(x - t(U) %*% c) #9
      }
    }
  } else { # Conditional sampling
    xM <- matrix(xM, ncol = 1)
    U <- setdiff(1:ncol(Q), M)
    QUU <- Q[U, U]
    QUM <- Q[U, M]
    if(missing(constr)) {
      if (use.inla) inla.qsample(n = n, Q = QUU, b = - QUM %*% (xM - mu[M]), mu = mu[U], logdens =  FALSE)
      else mu[U] + gmrf(n = n, Q = QUU, mu = t(-QUM %*% (xM - mu[M])))
    } else {
      if (use.inla) inla.qsample(n = n, Q = QUU, b = - QUM %*% (xM - mu[M]), mu = mu[U], constr = list(A = matrix(constr$A[,U], ncol = length(U)), e = constr$e), logdens =  FALSE)  
      else mu[U] + gmrf(n = n, Q = QUU, mu = t(-QUM %*% (xM - mu[M])), constr = list(A = rBind(constr$A[,U]), e = constr$e))
    }
  }
}
