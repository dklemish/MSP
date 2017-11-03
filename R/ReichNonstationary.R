library(stabledist)
library(evir)

set.seed(6)

n.sim <- 10000
alpha.center <- 0.9
alpha.surround <- 0.9

U <- rgev(n.sim, 1, alpha.center, alpha.center)
A <- matrix(rstable(4*n.sim, alpha.surround, 1, 1, 0, pm=0), ncol=4)

w.l <- 0.25
w.l.alpha <- w.l^(1/alpha.center)

sum.A <- apply(A*w.l.alpha, 1, sum)
X <- U*sum.A^alpha.center

hist(X, freq=FALSE, xlim=c(0,10), breaks=1000000)
lines(seq(0.01, 10, by=0.01), dgev(seq(0.01, 10, by=0.01), 1, 1, 1), col="red")