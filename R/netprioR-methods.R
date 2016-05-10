#' \code{\linkS4class{netprioR}} accessors
#' 
#' Retrieve network weights, prioritisation ranks and ROC curves
#'
#' @author Fabian Schmich
#' @rdname methods
#' @export
#' @param object A \code{\linkS4class{netprioR}} object 
#' @return Estimated network weights, prioritisation ranks or ROC curve
setGeneric(name = "weights", def = function(object) standardGeneric("weights"))
setGeneric(name = "ranks", def = function(object) standardGeneric("ranks"))
setGeneric(name = "ROC", def = function(object, ...) standardGeneric("ROC"))

#' @rdname methods
setMethod(f = "weights",
          signature = signature(object = "netprioR"),
          function(object) {
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

#' @rdname methods
#' @import dplyr
#' @inheritParams weights
#' @export
setMethod(f = "ranks",
          signature = signature(object = "netprioR"),
          function(object) {
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

#' @rdname methods
#' @importFrom pROC roc plot.roc
#' @inheritParams weights
#' @param true.labels True full set of underlying labels
#' @param plot Flag whether to plot the AUC curve
#' @export
setMethod(f = "ROC",
          signature = signature(object = "netprioR"),
          function(object, true.labels, plot = FALSE) {
            if(object@is.fitted) {
              stopifnot(length(object@model$Yimp) == length(true.labels))
              stopifnot((levels(object@labels)) == levels(true.labels))
              stopifnot(length(levels(true.labels)) == 2)
              unlabelled <- which(is.na(object@labels))
              ans <- roc(cases = object@model$Yimp[intersect(unlabelled, which(true.labels == levels(true.labels)[1]))], 
                         controls = object@model$Yimp[intersect(unlabelled, which(true.labels == levels(true.labels)[2]))])
              if (plot) plot.roc(ans, print.auc = TRUE, print.auc.x = 0.2, print.auc.y = 0.1)
              return(ans)
            } else {
              warning("No fitted model.")
            }
          }
)


#' Fit \code{\linkS4class{netprioR}} model
#'
#' @author Fabian Schmich
#' @export
#' @param object A \code{\linkS4class{netprioR}} object 
#' @return A \code{\linkS4class{netprioR}} object with fitted model
setGeneric(name = "fit", def = function(object, ...) standardGeneric("fit"))
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
