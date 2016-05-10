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

## Simulate data
dat <- simulate_labels(values = c(+1, -1), sizes = nmemb, nobs = c(nlabel/2, nlabel/2) %>% floor)
W <- list(TRU1 = simulate_network_scalefree(nmemb = nmemb, pclus = 0.8),
          TRU2 = simulate_network_scalefree(nmemb = nmemb, pclus = 0.8)
)
if (nrand > 0) {
  for (j in 1:nrand) {
    W[[sprintf("RAND%d", j)]] <- simulate_network_random(nmemb = nmemb, nnei = 1)  
  }          
}
W.norm <- W %>% lapply(., function(x) x / norm(x))
G <- lapply(W.norm, laplacian, norm = "none")
X <- simulate_phenotype(Y = dat$Y, meandiff = mdiff, sd = 1, pivot = 0.5)

## Fit netprioR
time.netprioR <- proc.time()
netprioR <- learn(Yobs = dat$Yobs, X = X, G = G, l = dat$l, u = dat$u, 
                  a = a, b = b, sigma2 = sigma2, tau2 = tau2, eps = eps, 
                  max.iter = max.iter, thresh = thresh, use.cg = use.cg, 
                  nrestarts = nrestarts, max.cores = max.cores, verbose = verbose)
time.netprioR <- (proc.time() - time.netprioR)

par(mfrow = c(2,1))
netprioR$logliks[] %>% plot(type = "l")
plot(netprioR$Yimp[dat$u])



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
