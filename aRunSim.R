remove(list = ls()) # clear environment
source("rfun.R"); source("rsim.R"); source("rlab_x.R");


# PARAMS ------------------------------------------------------------------


p <- 4
N <- 1e3
Matrix_0 <- matrix(rep(0, p^2), ncol = p)

Omega    <- matrix(runif(p^2, 1, 10), ncol = p)
Epsilon  <- 1e-2 * matrix(runif(p^2, 1, 10), ncol = p)


# sources and mixtures ----------------------------------------------------


z <- sim_good_sources(N, p)
x <- tvmix(z, Omega, Epsilon)

xc <- scale(x, scale = F)
res <- ltvsobi(xc)
plot.ts(res$S, ylim = c(-10, 10))
plot.ts(z)
plot(res$S[,4], ylim = c(-10,10))

summary(res$S) 

for (i in 1:p) {
  restored <- res$S[,i]
  plot(restored, ylim = quantile(restored, probs = c(0.1, .9)))
}


# fig_mixing <- list(); fig_mixing$source <- z
# fig_mixing$tvmix<- tvmix(z, Omega, Epsilon); fig_mixing$mix <- tvmix(z, Omega, Matrix_0); fig_mixing$unmix<- JADE::SOBI(fig_mixing$mix)$S
# save(fig_mixing, file = paste0(getwd(), "/thesis/thesis.rdata"))

x <- scale(x, center = T)

res0 <- JADE::SOBI(x)
res1 <- ltvsobi(x)
res2 <- ltvsobi2(x)
res9 <- tvsobi(x)
SIR_all(res0, Omega, Epsilon, z)
SIR_all(res1, Omega, Epsilon, z)
SIR_all(res2, Omega, Epsilon, z)
SIR_all(res9, Omega, Epsilon, z)

