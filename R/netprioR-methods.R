#' Retrieve network weights and prioritisation ranks from \code{\linkS4class{netprioR}} objects
#'
#' @author Fabian Schmich
#' @rdname methods
#' @export
#' @param object A \code{\linkS4class{netprioR}} object 
#' @return Estimated network weights or prioritisation ranks
#' @seealso \code{\link{ranks}}
#' @seealso \code{\link{weights}}
setGeneric(name = "weights", def = function(object) standardGeneric("weights"))

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


setGeneric(name = "ranks", def = function(object) standardGeneric("ranks"))
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
