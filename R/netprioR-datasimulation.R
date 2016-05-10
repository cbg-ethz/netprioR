#' Simulate scalefree networks
#'
#' Simulate scale free networks for predefined number of members for each of
#' two groups and a parameter pclus that determines how strictly distinct the groups
#' are
#' 
#' @author Fabian Schmich
#' @import Matrix
#' @export 
#' @param nmemb Vector of numbers of members per group
#' @param pclus Scalar in [0, 1] determining how strictly distinct groups are
#' @return Adjacency matrix
simulate_network_scalefree <- function(nmemb = c(100, 100), pclus = 1) {
  N <- sum(nmemb)
  names <- paste(rep(LETTERS[1:length(nmemb)], nmemb), sapply(nmemb, function(x) 1:x), sep = "")
  X <- Matrix(0, nrow = N, ncol = N, dimnames = list(names, names))  
  for (r in 1:length(nmemb)) {
    offs <- ifelse(r == 1, 0, nmemb[1:(r-1)] %>% sum)
    grpset <- seq(offs + 1, offs + nmemb[r])
    for (i in 1:length(grpset)) {
      if (runif(1) <= pclus) {
        gset <- setdiff(grpset, offs + i) # all nodes in the group but itself
      } else {
        gset <- setdiff(1:N, grpset) # all nodes in the other groups
      }
      if (all(colSums(X[gset,gset]) == 0)) { # 1st vertex
        at <- sample(gset, size = 1, replace = FALSE)
      } else { # preferential attachment
        at <- sample(gset, 
                     size = 1,
                     prob = colSums(X[gset, gset]) / sum(colSums(X[gset, gset])),
                     replace = FALSE)
      }
      X[at, offs + i] <- X[offs + i, at] <- 1
    }
  }
  if(any(colSums(X) == 0)) stop("Created un-attached vertex")
  return(X)
}

#' Simulate random networks with predefined number of members for each
#' of the two groups and the number of neighbours for each node
#' 
#' @author Fabian Schmich
#' @import Matrix
#' @export
#' @param nmemb Vector of number of members for each group
#' @param nnei Number of neighbours for each node
#' @return Adjacency matrix of graph
simulate_network_random <- function(nmemb = c(100, 100), nnei = 1) {
  N <- sum(nmemb)
  names <- paste(rep(LETTERS[1:length(nmemb)], nmemb), sapply(nmemb, function(x) 1:x), sep = "")
  X <- Matrix(0, nrow = N, ncol = N, dimnames = list(names, names))  
  for (i in 1:nrow(X)) {
    neis <- sample(setdiff(1:N, i), size = nnei, replace = FALSE)
    X[i,neis] <- X[neis,i] <- 1
  }
  return(X)
}


#' Simulate labels
#' 
#' @author Fabian Schmich
#' @import dplyr
#' @export
#' @param values Vector of labels for groups
#' @param sizes Vector of group sizes
#' @param nobs Vector of number of observed labels per group
#' @return List of Y, Yobs and indices for labeled instances
simulate_labels <- function(values, sizes, nobs) {
  stopifnot(length(sizes) == length(values) & length(nobs) == length(values))
  Y <- sapply(1:length(values), function(i) rep(values[i], sizes[i])) %>% as.numeric
  l <- sapply(1:length(values), function(i) {
    if (i == 1) {
      sampfrom <- (1:sizes[i])
    } else {
      sampfrom <- (sizes[1:(i-1)] %>% sum + 1):(sizes[1:(i-1)] %>% sum + sizes[i])
    }
    sample(sampfrom, nobs[i]) 
  }) %>% unlist %>% as.numeric %>% sort
  u <- setdiff(1:length(Y), l)
  #   Y <- c(Y[l], Y[u])
  #   l <- 1:sum(nobs)
  #   u <- setdiff(1:length(Y), l)
  Yobs <- Y
  Yobs[u] <- 0
  return(list(Y = Y, Yobs = Yobs, N = length(Y), l = l, u = u))
}

#' Simulate phenotypes correlated to labels pivoted into two groups
#' 
#' @author Fabian Schmich
#' @import Matrix
#' @export
#' @param Y Vector of labels
#' @param meandiff difference of means between positive and negative groups
#' @param sd Standard deviation of the phenotype
#' @param pivot Pivot point to split labels into groups
#' @return Simulated phenotype
simulate_phenotype <- function(Y, meandiff, sd, pivot = 0) {
  X <-  rep(NA, length(Y)) %>% cbind
  s0 <- which(Y < pivot)
  b0 <- which(Y >= pivot)
  X[s0] <- rnorm(length(s0), -meandiff/2, sd)
  X[b0] <- rnorm(length(b0), +meandiff/2, sd)
  return(Matrix(X))
}
