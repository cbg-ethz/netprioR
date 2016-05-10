#' @export
if (!isGeneric("summary")) {
  setGeneric(name = "summary",
             def = function(object, ...) {
               standardGeneric("summary")
             },
             package = "netprioR"
  )  
}
setMethod(f = "summary",
          signature = signature(object = "netprioR"),
          definition = function(object) {
            show(object)
            cat("Model: ")
            if (object@is.fitted) {
              cat("\nLikelihood[log]:", tail(object@model$logliks, n = 1), "\n")
              cat("Fixed effects:", object@model$beta, "\n")
              cat("Network weights:\n")
              print(data.frame(Network = names(object@model$W), Weight = object@model$W), 
                    quote = FALSE, row.names = FALSE)
            } else {
              cat("not fitted.\n")
            }    
          }
)

setMethod(f = "show",
          signature = signature(object = "netprioR"),
          definition = function(object) {
            cat("#Genes:", length(object@labels), "\n")
            cat("#Networks: ", length(object@networks), "\n")
            cat("#Phenotypes:", ncol(object@phenotypes), "\n")
            cat("#Labels:", length(which(!is.na(object@labels))), "\n")
            cat("Classes:", levels(object@labels), "\n")
          }
)
