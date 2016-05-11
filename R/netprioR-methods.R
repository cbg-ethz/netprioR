#' Retrieve network weights
#'
#' @author Fabian Schmich
#' @rdname weights-methods
#' @export
#' @param object A \code{\linkS4class{netprioR}} object 
#' @param ... Additional arguments
#' @return Estimated network weights
setGeneric(name = "weights", def = function(object) standardGeneric("weights"))
#' @rdname weights-methods
#' @examples 
#' data(simulation)
#' weights(simulation$model)
setMethod(f = "weights",
          signature = signature(object = "netprioR"),
          function(object) {
            Network <- Weight <- NULL
            if(object@is.fitted) {
              ans <- data.frame(Network = names(object@model$W), Weight = object@model$W) %>%
                tbl_df %>%
                arrange(desc(Weight))
              rownames(ans) <- c()
              return(ans)
            } else {
              warning("No fitted model.")
            }
          }
)

#' Retrieve ranked prioritisation list
#' 
#' @author Fabian Schmich
#' @rdname ranks-methods
#' @import dplyr
#' @param object A \code{\linkS4class{netprioR}} object 
#' @return Ranked list of prioritised genes
#' @export
setGeneric(name = "ranks", def = function(object) standardGeneric("ranks"))
#' @rdname ranks-methods
#' @examples 
#' data(simulation)
#' ranks(simulation$model)
setMethod(f = "ranks",
          signature = signature(object = "netprioR"),
          function(object) {
            Score <- Rank <- Id <- NULL
            if(object@is.fitted) {
              unlabelled <- which(is.na(object@labels))
              # Try to find labels somewhere
              if (!is.null(rownames(object@phenotypes))) {
                ids <- rownames(object@phenotypes)
              } else if (!is.null(names(object@labels))) {
                ids <- names(object@labels)
              } else {
                ids <- 1:length(object@labels)  
              }
              data.frame(Id = ids[unlabelled], 
                         Score = object@model$Yimp[unlabelled]) %>%
                tbl_df %>%
                arrange(desc(Score)) %>%
                mutate(Rank = 1:length(unlabelled)) %>%
                select(Rank, Id, Score)
            } else {
              warning("No fitted model.")
            }
          }
)

#' Compute ROC curve from netprioR model and true labels
#'
#' @author Fabian Schmich
#' @rdname ROC-methods
#' @importFrom pROC roc plot.roc
#' @param object A \code{\linkS4class{netprioR}} object 
#' @param true.labels True full set of underlying labels
#' @param plot Flag whether to plot the AUC curve
#' @param ... Additional arguments
#' @return ROC curve with AUC
#' @export
setGeneric(name = "ROC", def = function(object, ...) standardGeneric("ROC"))
#' @rdname ROC-methods
#' @examples 
#' data(simulation)
#' ROC(simulation$model, true.labels = simulation$labels.true)
setMethod(f = "ROC",
          signature = signature(object = "netprioR"),
          function(object, true.labels, plot = FALSE, ...) {
            if(object@is.fitted) {
              stopifnot(length(object@model$Yimp) == length(true.labels))
              stopifnot((levels(object@labels)) == levels(true.labels))
              stopifnot(length(levels(true.labels)) == 2)
              unlabelled <- which(is.na(object@labels))
              roccurve(labels.pred = object@model$Yimp[unlabelled], labels.true = true.labels[unlabelled], plot = plot, ...)
              # ans <- roc(cases = object@model$Yimp[intersect(unlabelled, which(true.labels == levels(true.labels)[1]))], 
              #            controls = object@model$Yimp[intersect(unlabelled, which(true.labels == levels(true.labels)[2]))],
              #            direction = ">")
              # if (plot) plot.roc(ans, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.1, ...)
              # return(ans)
            } else {
              warning("No fitted model.")
            }
          }
)


#' Fit \code{\linkS4class{netprioR}} model
#'
#' @author Fabian Schmich
#' @rdname fit-methods
#' @export
#' @param object A \code{\linkS4class{netprioR}} object 
#' @param refit Flag whether to overwrite existing fit
#' @param ... Additional arguments
#' @return A \code{\linkS4class{netprioR}} object with fitted model
setGeneric(name = "fit", def = function(object, ...) standardGeneric("fit"))
#' @rdname fit-methods
#' @examples 
#' data(simulation)
#' np <- netprioR(networks = simulation$networks,
#'                phenotypes = simulation$phenotypes,
#'                labels = simulation$labels.obs,
#'                model.fit = FALSE)
#' summary(np)
#' np <- fit(np, nrestarts = 1, verbose = FALSE)
#' summary(np)
setMethod(f = "fit",
          signature = signature(object = "netprioR"),
          function(object, refit = FALSE, ...) {
            if(!object@is.fitted) {
              netprioR(networks = object@networks,
                       phenotypes = object@phenotypes,
                       labels = object@labels,
                       fit.model = TRUE,
                       ...)
            } else if (refit) {
              netprioR(networks = object@networks,
                       phenotypes = object@phenotypes,
                       labels = object@labels,
                       fit.model = TRUE,
                       ...)
            } else {
              warning("Set refit = TRUE, if existing fit should be overwritten.")
              return(object)
            }
          }
)


#' Retrieve elements of a fitted \code{\linkS4class{netprioR}} object
#'
#' @author Fabian Schmich
#' @rdname getter-methods
#' @export
#' @param object A \code{\linkS4class{netprioR}} object 
#' @return Model of a fitted \code{\linkS4class{netprioR}} object
setGeneric(name = "model", def = function(object) standardGeneric("model"))
#' @rdname getter-methods
setMethod(f = "model",
          signature = signature(object = "netprioR"),
          function(object) {
            if(object@is.fitted) {
              return(object@model)
            } else {
              warning("No fitted model.")
              return(object)
            }
          }
)

#' @rdname getter-methods
#' @inheritParams model
#' @export
#' @return Phenotypes of a fitted \code{\linkS4class{netprioR}} object
setGeneric(name = "phenotypes", def = function(object) standardGeneric("phenotypes"))
#' @rdname getter-methods
setMethod(f = "phenotypes",
          signature = signature(object = "netprioR"),
          function(object) {
            return(object@phenotypes)
          }
)

#' @rdname getter-methods
#' @inheritParams model
#' @export
#' @return Networks of a fitted \code{\linkS4class{netprioR}} object
setGeneric(name = "networks", def = function(object) standardGeneric("networks"))
#' @rdname getter-methods
setMethod(f = "networks",
          signature = signature(object = "netprioR"),
          function(object) {
            return(object@networks)
          }
)
