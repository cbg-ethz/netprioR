#' netprioR
#' 
#' Class that represents a netprioR model.
#' 
#' @author Fabian Schmich
#' @name netprioR-class
#' @rdname netprioR-class
#' @aliases netprioR
#' 
#' @import dplyr
#' @exportClass netprioR
#'
#' @slot networks List of NxN adjacency matrices of gene-gene similarities
#' @slot phenotypes Matrix of dimension NxP containing covariates
#' @slot labels Vector of Nx1 labels for all genes. NA if no label available. 
#' @slot is.fitted Flag indicating if model is fitted
#' @slot model List containing estimated parameters and imputed missing data
#' 
setClass(Class = "netprioR",
         representation = representation(
           networks = "list",
           phenotypes = "matrix",
           labels = "factor",
           is.fitted = "logical",
           model = "list"
         ),
         validity=function(object) {
           return(TRUE)
           # TODO: implement validity check
         }
)

#' @rdname netprioR-class
#' @exportMethod netprioR
#'
#' @param networks List of NxN adjacency matrices of gene-gene similarities
#' @param phenotypes Matrix of dimension NxP containing covariates
#' @param labels Vector of Nx1 labels for all genes (NA if no label available)
#' @param fit.model Indicator whether to fit the model
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
#' @param ... Additional arguments
#' @return A \code{\linkS4class{netprioR}} object
setGeneric("netprioR",
           function(networks, phenotypes, labels, ...) {
             standardGeneric("netprioR")
           })
#' @rdname netprioR-class
#' @examples 
#' \donttest{ # runs long-ish
#' data(simulation)
#' np <- netprioR(networks = simulation$networks,
#'                phenotypes = simulation$phenotypes,
#'                labels = simulation$labels.obs,
#'                fit.model = TRUE)
#' summary(np)          
#' }
setMethod("netprioR",
          signature = signature(networks = "list", 
                                phenotypes = "matrix",
                                labels = "factor"),
          function(networks, 
                   phenotypes, 
                   labels,
                   fit.model = FALSE,
                   a = 0.1,
                   b = 0.1,
                   sigma2 = 0.1,
                   tau2 = 100,
                   eps = 1e-10,
                   max.iter = 500,
                   thresh = 1e-6,
                   use.cg = FALSE,
                   thresh.cg = 1e-6,
                   nrestarts = 5,
                   max.cores = detectCores(),
                   verbose = TRUE,
                   ...){
            stopifnot(length(levels(labels)) == 2)
            if (fit.model) {
              labelled <- which(!is.na(labels))
              unlabelled <- setdiff(1:length(labels), labelled)
              Yobs <- factor(labels, labels = c(-1, +1)) %>% as.character %>% as.numeric
              G <- lapply(networks, function(x) {
                x <- x / norm(x)
                return(laplacian(x, norm = "none"))
              })
              model <- learn(Yobs = Yobs,
                             G = G,
                             X = phenotypes,
                             l = labelled,
                             u = unlabelled,
                             a = a,
                             b = b,
                             sigma2 = sigma2,
                             tau2 = tau2,
                             eps = eps,
                             max.iter = max.iter,
                             thresh = thresh,
                             use.cg = use.cg,
                             nrestarts = nrestarts,
                             max.cores = max.cores,
                             verbose = verbose)
              
              new("netprioR",
                  networks = networks,
                  phenotypes = phenotypes,
                  labels = labels,
                  model = model,
                  is.fitted = TRUE)
            } else {
              new("netprioR",
                  networks = networks,
                  phenotypes = phenotypes,
                  labels = labels,
                  model = list(),
                  is.fitted = FALSE)
            }
          }
)

#' Plot method for \code{\linkS4class{netprioR}} objects
#' 
#' @author Fabian Schmich
#' @import ggplot2
#' @import dplyr
#' @importFrom gridExtra grid.arrange
#' @export
#' @method plot netprioR
#' 
#' @param x A \code{\linkS4class{netprioR}} object
#' @param which Flag for which plot should be shown, options: weights, lik, scores, all
#' @param ... Additional paramters for plot
#' @return Plot of the weights, likelihood, ranks, or all three
#' @examples 
#' data(simulation)
#' plot(simulation$model)
plot.netprioR <- function(x, which = c("all", "weights", "lik", "scores"), ...) {
  which <- match.arg(which)
  Weight <- Network <- Iteration <- Loglik <- Score <- Rank <- Id <- NULL
  if (x@is.fitted) {
    pl.weights <- ggplot(weights(x) %>% 
                           mutate(Weight = Weight / sum(Weight)) %>%
                           arrange(desc(Weight)) %>%
                           mutate(Network = factor(Network, levels = Network)), 
                         aes(x = Network, y = Weight)) + 
      geom_bar(stat = "identity") + 
      ylab("Relative weight") + 
      ggtitle("Network weights") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    pl.lik <-ggplot(data.frame(Iteration = 2:length(x@model$logliks), Loglik = x@model$logliks[-1]), aes(x = Iteration, y = Loglik)) + 
      geom_line() + 
      scale_y_continuous(labels = function(x) format(x, nsmall = 2, scientific = TRUE)) +
      ylab("Likelihood [log]") + 
      ggtitle("EM iterations") +
      theme_bw()
    pl.scores <- ggplot(ranks(x), aes(x = Score)) + geom_histogram(binwidth = 0.1) + 
      ylab("Count") +
      ggtitle("Prioritisation") +
      theme_bw()
    switch(which,
           "weights" = pl.weights,
           "lik" = pl.lik,
           "scores" = pl.scores,
           "all" = grid.arrange(pl.lik, pl.weights, pl.scores, ncol = 1)
    )
  } else {
    warning("No model fitted.")
  }
}
