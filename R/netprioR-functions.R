#' Fit netprioR model
#' 
#' Infer parameters and hidden data using the EM algorithm of netprioR
#' 
#' @author Fabian Schmich
#' @import Matrix
#' @import dplyr
#' @importFrom stats runif rnorm
#' @importFrom sparseMVN dmvn.sparse rmvn.sparse
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores
#' @importFrom foreach foreach %dopar%
#' @param Yobs Observed labels (NA, if not observed)
#' @param X Phenotypes
#' @param G Graph Laplacians
#' @param l Indices of labelled instances
#' @param u Indices of unlabelled instances
#' @param a Shape parameter of Gamma prior for W
#' @param b Scale parameter of Gamma prior for W
#' @param sigma2 Cariance for Gaussian labels
#' @param tau2 Variance for Gaussian prior for beta
#' @param eps Small value added to diagonal of Q in order to make it non-singular
#' @param max.iter Maximum number of iterations for EM
#' @param thresh Threshold for termination of EM with respect to change in parameters
#' @param use.cg Flag whether to use conjugate gradient instead of exact computation of expectations
#' @param thresh.cg Threshold for the termination of the conjugate gradient solver
#' @param nrestarts Number of restarts for EM
#' @param max.cores Maximum number of cores to use for parallel computation
#' @param verbose Print verbose output
#' @return List containing: Predicted labels Yhat and inferred parameters W and beta
learn <- function(Yobs, X, G, l, u, a = 0.1, b = 0.1, sigma2 = 1, tau2 = 10, eps = 1e-11, max.iter = 500, 
                  thresh = 1e-3, use.cg = TRUE, thresh.cg = 1e-5, nrestarts =  5, max.cores = detectCores(), verbose = FALSE) {
  ## Initialise variables and settings
#   options(digits = 17)
  K <- length(G)
  N <- length(Yobs)
  P <- ncol(X)
  R.new <- rep(0, N)
  Yimp <- vector(length = N)
  loglik <- loglik.new <- -.Machine$double.xmin
  conv <- 1
  Smat <- sigma2 * Diagonal(N)
  Tmat <- tau2 * Diagonal(P)
  logliks <- c()
  
  ## Initialise multicore
  if(nrestarts > 1 & max.cores > 1) {
    registerDoParallel(cores = min(nrestarts, max.cores))
  }
  
  EM.runs <- foreach (runs = 1:nrestarts) %dopar% {
    W <- W.new <- runif(K, min = 1, max = 1000) #rgamma(n = K, shape = a, rate = b) #rep(1, K)
    beta <- beta.new <- rmvn.sparse(n = P, mu = rep(0, P), CH = Cholesky(Tmat), prec = FALSE)
    Q <- sapply(1:K, function(k) W[k] * G[[k]]) %>% Reduce("+", .) + eps * Diagonal(N)
    for (i in 1:max.iter) {
      ## Expectation Step
      ## E[R]
      H <- Yobs - X %*% beta
      if (use.cg) {
        E_RL <- solve(solve(Smat[l,l]) + Q[l,l]) %*% Q[l,l] %*% H[l]
        E_RU <- conjugate_gradient(A = Q[u,u], b = -Q[u,l] %*% E_RL, threshold = thresh.cg, verbose = FALSE)
      } else {
        E_RL <- solve(solve(Smat[l,l]) + Q[l,l]) %*% Q[l,l] %*% H[l]
        E_RU <- solve(Q[u,u]) %*% -Q[u,l] %*% E_RL
      }   
      R.new[u] <- E_RU %>% as.numeric
      R.new[l] <- E_RL %>% as.numeric
      
      ## E[Yu]
      E_YU = (R.new + X %*% beta)[u] %>% as.numeric
      Yimp[l] <- Yobs[l]
      Yimp[u] <- E_YU
      
      ## Maximisation step: W, beta
      W.new <- sapply(1:K, function(k) ((a - 1 + N/2) / (b + 0.5 * (R.new %*% G[[k]] %*% R.new))) %>% as.numeric)
      beta.new <- solve(t(X[l, drop = FALSE]) %*% solve(Smat[l,l]) %*% X[l, drop = FALSE] + solve(Tmat)) %*% t(X[l, drop = FALSE]) %*% solve(Smat[l,l]) %*% (Yimp[l] - R.new[l])
      #       beta.new <- solve(t(X) %*% solve(Smat) %*% X + solve(Tmat)) %*% t(X) %*% solve(Smat) %*% (Yimp - R.new)
      
      ## Compute precision matrix
      Q <- sapply(1:K, function(k) W.new[k] * G[[k]]) %>% Reduce("+", .) + eps * Diagonal(N)
      
      ## Compute log liklihood up to constants
      loglik.new <- c(- t(Yimp - (X %*% beta.new + R.new))[l] %*%
                        solve(Smat[l,l]) %*%
                        (Yimp - (X %*% beta.new + R.new))[l],
                      - t(R.new) %*% Q %*% R.new,
                      - t(beta.new) %*% solve(Tmat) %*% beta.new,
                      (a - 1) * log(W.new) - b * W.new) %>% sapply(as.numeric) %>% sum
      
      # loglik.new <- c(dmvn.sparse(x = Yimp[l], mu = (X %*% beta.new + R.new)[l] %>% as.numeric, CH = Cholesky(Smat[l,l]), prec = FALSE),
      #                 dmvn.sparse(x = R.new, mu = rep(0, N), CH = Cholesky(Q), prec = TRUE),
      #                 dmvn.sparse(x = beta.new %>% as.numeric, mu = rep(0, P), CH = Cholesky(Tmat), prec = FALSE),
      #                 dgamma(x = W.new, shape = a, rate = b, log = TRUE)) %>% sum

      logliks <- c(logliks, loglik.new)
      
      ## Check convergence
      conv <- abs((loglik.new - loglik) / loglik)
      if (verbose) {
        cat(i, "iteration: ", conv, "\n")
        #cat("W: ", W.new %>% as.numeric, "\n")
        cat("beta: ", beta.new %>% as.numeric, "\n")
        print(data.frame(G = names(G), W = W.new) %>% tbl_df %>% arrange(desc(W)) %>% data.frame)
        cat("loglik: ", loglik.new %>% as.numeric, "\n\n")
      }
      
      ## Terminate if loglik doesn't change more than given eps
      if(conv < thresh) {
        break
      } else {
        ## Update parameters
        W <- W.new
        beta <- beta.new
        loglik <- loglik.new
      }
      
    }
    W.new <- as.numeric(W.new)
    names(W.new) <- names(G)
    beta.new <- as.numeric(beta.new)
    names(beta.new) <- colnames(X)
    return(list(Yimp = Yimp,
                R = R.new %>% as.numeric,
                W = W.new,
                beta = beta.new,
                logliks = logliks))
  }
  ## Select run with maximum liklihood
  maxlikind <- lapply(EM.runs, function(x) {
    tail(x$logliks, n = 1)
  }) %>% unlist %>% which.max
  return(EM.runs[[maxlikind]])
}

#' Graph Laplacian
#' 
#' Compute the Laplacian matrix of a graph given its adjacency matrix
#'
#' @author Fabian Schmich
#' @import Matrix
#' @param x Adjacency matrix
#' @param norm Type of normalisation
#' @return Laplacian matrix
laplacian <- function(x, norm = c("none", "sym", "asym")) {
  norm <- match.arg(norm)
  stopifnot(nrow(x) == ncol(x)) 
  d <- rowSums(x)
  dsqrt <- sapply(d, function(x) ifelse(x == 0, 0, 1/sqrt(x)))
  D <- Diagonal(x = d)
  if (!is.null(dimnames(x))) dimnames(D) <- dimnames(x)
  L <- D - x
  switch(norm,
         none = L,
         sym = Diagonal(x = dsqrt) %*% L %*% Diagonal(x = dsqrt),
         asym = solve(D) %*% L)
}


#' Conjugate Gradient Solver
#' 
#' Solves linear equation systems iteratively
#' 
#' @author Fabian Schmich
#' @import Matrix
#' @param A Matrix
#' @param b Coefficients
#' @param x0 Starting solution
#' @param threshold Termination threshold
#' @param verbose Show iterative progress
#' @return Solution for equation system
conjugate_gradient <- function(A, b, x0 = rep(0, ncol(A)), threshold = 1e-15, verbose = FALSE) {
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

#' Class Mass Normalization (CMN) from Zhu et al., 2003
#' 
#' @author Fabian Schmich
#' @import Matrix
#' @importFrom stats pnorm
#' @param yhat Response for labeled (l) and unlabeld (u) genes
#' @param l Indices of labeled genes
#' @param u Indices of unlabeled genes
#' @return Class normalized yhat
cmn <- function(yhat, l, u) {
  fL <- cbind(pnorm(yhat[l]), 1 - pnorm(yhat[l]))
  fU <- cbind(pnorm(yhat[u]), 1 - pnorm(yhat[u]))
  N <- nrow(fU) + nrow(fL)
  q <- (colSums(fL) + 1)
  print(q)
  ans <- fU * Matrix(rep(q / colSums(fU), each = N - nrow(fL)), nrow = N - nrow(fL))
  apply(ans, 1, function(x) ifelse(x[1] > x[2], 1, 0))
}

#' bandwidth
#' 
#' Compute the bandwidth of a matrix
#'
#' @author Fabian Schmich
#' @import Matrix
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
#' @import Matrix
#' @param x Input matrix
#' @return Band matrix
cuthill_mckee <- function(x) {
  degs <- data.frame(Idx=1:ncol(x), NonZero=apply(x, 1, function(x) length(which(x != 0))))
  R <- degs$Idx[which.min(degs$NonZero)]
  i <- 1
  for (i in 1:ncol(x)) {
    Ai <- setdiff(which(x[R[i],] != 0), R)
    if (length(Ai) > 0) {
      Ai <- Ai[order(degs$NonZero[Ai], decreasing = FALSE)]
      R <- append(R, Ai)
    } else {
      R <- append(R, degs$Idx[-R][which.min(degs$NonZero[-R])])
    }
    i <- i + 1
  }
  rR <- rev(R)
  return(x[rR, rR])
}

#' Normalise kernel
#' 
#' adopted from GeneMania, Mostafavi et al, 2009
#' 
#' @author Fabian Schmich
#' @import Matrix
#' @param x kernel
#' @return Normalised kernel
norm_kern <- function(x) {
  x <- x - diag(x)
  nz.ind <- which(x != 0, arr.ind = TRUE)
  nz.val <- x[nz.ind]
  colsums.sqrt <- sapply(colSums(x), function(i) 1/sqrt(i))
  for (i in 1:nrow(nz.ind)) {
    rr <- colsums.sqrt[nz.ind[i,1]]
    cc <- colsums.sqrt[nz.ind[i,2]]
    nz.val[i] <- nz.val[i] * rr * cc
  }
  x[nz.ind] <- nz.val
  return(x)
}
