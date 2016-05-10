require(dplyr)
require(Matrix)

# Settings
N <- 1000
nmemb <-  c(N/2, N/2) %>% floor
nlabel <- 50
nrand <- 1
mdiff <- 0.5
a <- b <- 0.01
sigma2 <- 0.1
tau2 <- 100
eps <- 1
max.iter <- 1000
thresh <- 1e-6
use.cg <- FALSE
nrestarts <- 1
max.cores <- 1
verbose <- TRUE

networks <- W
phenotypes <- as.matrix(X)
true.labels <- dat$Y %>% sapply(function(x) ifelse(x == 0, NA, x)) %>% factor(levels = c(-1, 1), labels = c("Negative", "Positive"))
labels <- dat$Yobs %>% sapply(function(x) ifelse(x == 0, NA, x)) %>% factor(levels = c(-1, 1), labels = c("Negative", "Positive"))
np <- netprioR(networks = networks, phenotypes = phenotypes, labels = labels, nrestarts = 1, thresh = 1e-6, a = 0.1, b = 0.1, fit.model = TRUE)
summary(np)
np

plot(np)
plot(np, which = "lik")
plot(np, which = "weight")
plot(np, which = "scores")

ROC(np, true.labels, plot = TRUE)


np2 <- netprioR(networks = networks, phenotypes = phenotypes, labels = labels, fit.model = FALSE)
np2
summary(np2)
plot(np2)
np2 <- fit(np2, verbose = TRUE, thresh = 1e-3)
plot(np2)
np2 <- fit(np2, thresh = 1e-4, refit = TRUE)
plot(np2)
