
N <- 10
Q <- Matrix(rbinom(N^2, 1, 0.2), nrow = N)
Q <- forceSymmetric(Q)
Q <- laplacian(Q)
all(eigen(Q)$values > 0)
chol(Q)



g = system.file("demodata/germany.graph", package="INLA")
G = inla.graph2matrix(g)
lp <- runif(n = ncol(G))
Q <- pmat(genes = 1:ncol(G), interactions = G, label.priors = list("P" = lp, "N" = 1 - lp), theta = 1)

M <- 1:20
xM <- c(rep(-5, 10), rep(5, 10))
A <- matrix(rep(1, ncol(Q)), ncol = ncol(Q))
constr = list(A = A, e = 0)
x <- gmrf(n = 1000, Q = Q, mu = 0, constr = list(A = A, e = 0), M = M, xM = xM, use.inla = FALSE)


g = system.file("demodata/germany.graph", package="INLA")
Q = inla.graph2matrix(g)
diag(Q) <- dim(Q)[1]

A <- matrix(rep(1, ncol(Q)), ncol = ncol(Q))
constr = list(A = A, e = 0)

gmrf(n = 1000, Q = Q, mu = 0, constr = list(A = A, e = 0), M = 1:2, xM = c(-5,5), use.inla = FALSE)
gmrf(n = 10, Q = Q, mu = 0, M = 1:2, xM = c(-5,5), use.inla = TRUE)
gmrf(n = 10, Q = Q, mu = 0, use.inla = TRUE)



M <- 1:20
xM <- c(rep(-5, 10), rep(5, 10))
a <- gmrf(n = 1000, Q = Q, mu = 0, constr = list(A = A, e = 0), M = M, xM = xM, use.inla = FALSE)
b <- gmrf(n = 1000, Q = Q, mu = 0, constr = list(A = A, e = 0), M = M, xM = xM, use.inla = TRUE)


x <- apply(a, 1, mean)
y <- apply(b, 1, mean)

plot(x, y)
