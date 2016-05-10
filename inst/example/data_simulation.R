## Load libraries
require(dplyr)

## Settings
set.seed(17)
N <- 1000 # number of genes
nlabel <- 50 # nubmer of labelled genes

## Simulate data
networks <- list(LOW_NOISE1 = simulate_network_scalefree(nmemb = c(N/2, N/2), pclus = 0.8),
                 LOW_NOISE2 = simulate_network_scalefree(nmemb = c(N/2, N/2), pclus = 0.8),
                 HIGH_NOISE = simulate_network_random(nmemb = c(N/2, N/2), nnei = 1)
)
labels <- simulate_labels(values = c("Positive", "Negative"), sizes = c(N/2, N/2), nobs = c(nlabel/2, nlabel/2)) 
phenotypes <- simulate_phenotype(labels.true = labels$labels.true, meandiff = 0.5, sd = 1)
np <- netprioR(networks = networks, phenotypes = phenotypes, labels = labels$labels.obs, fit.model = TRUE)
simulation <- c(labels[c("labels.true", "labels.obs")], networks = list(networks), phenotypes = list(phenotypes), model = np)
save(simulation, file = "data/simulation.rdata")
