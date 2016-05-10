#' Package: netprioR
#' 
#' This package provides a model for semi-supervised prioritisation of genes integrating network data, 
#' phenotypes and additional prior knowledge about TP and TN gene labels.
#' 
#' @name netprioR-package
#' @author Fabian Schmich | Computational Biology Group, ETH Zurich | \email{fabian.schmich@@bsse.ethz.ch}
#' @docType package
#' @keywords package
#' @references Fabian Schmich et. al (2016).
#' @import methods
#' @importFrom graphics plot
globalVariables(".")

#' Example data: Simulated networks, phenotypes and labels for N = 1000 genes
#' 
#' The data set contains simulated data for N = 1000 genes and P = 1 (univariate) phenotypes. 
#' The list of networks contains 2 low noise networks and two high noise networks. The class
#' labels are "Positive" and "Negative".
#' 
#' The code used to simluate the data can be found in system.file("example", "data_simulation.R", package = "netprioR")
#' 
#' @docType data
#' @name simulation
#' @usage data(simulation)
NA
