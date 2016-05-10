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
#' @param labels Vector of Nx1 labels for all genes. NA if no label available. 
#' @param ... Additional arguments
#' @return A \code{\linkS4class{netprioR}} object
setGeneric("netprioR",
           function(networks, phenotypes, labels, ...) {
             standardGeneric("netprioR")
           })
#' @rdname netprioR-class
setMethod("netprioR",
          signature = signature(networks = "list", 
                                phenotypes = "matrix",
                                labels = "factor"),
          function(networks, 
                   phenotypes, 
                   labels,
                   a = 0.1,
                   b = 0.1,
                   sigma2 = 1,
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
          }
)

#' Plot method for \code{\linkS4class{netprioR}} objects
#' 
#' @author Fabian Schmich
#' @import ggplot2
#' @import dplyr
#' @export
#' @method plot netprioR
#' 
#' @param x A \code{\linkS4class{netprioR}} object
#' @param which Flag for which plot should be shown, options: weights, lik, scores
#' @param ... Additional paramters for plot
#' @return Plot of the weights, likelihood, or ranks
plot.netprioR <- function(x, which = c("weights", "lik", "scores"), ...) {
  which <- match.arg(which)
  if (x@is.fitted) {
    switch(which,
           "weights" = ggplot(weights(x) %>% mutate(Weight = Weight / sum(Weight)), aes(x = Network, y = Weight)) + 
             geom_bar(stat = "identity") + 
             ylab("Relative weight") + 
             coord_flip() + 
             theme_bw(),
           "lik" = ggplot(data.frame(Iteration = 2:length(x@model$logliks), Loglik = x@model$logliks[-1]), aes(x = Iteration, y = Loglik)) + 
             geom_line() + 
             scale_y_continuous(labels = function(x) format(x, nsmall = 2, scientific = TRUE)) +
             ylab("Likelihood [log]") + 
             theme_bw(),
           "scores" = ggplot(ranks(x), aes(x = Score)) + geom_histogram(binwidth = 0.1) + 
             ylab("Count") +
             theme_bw()
    )
  } else {
    warning("No model fitted.")
  }
}
