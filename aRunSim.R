remove(list = ls()) # clear environment
source("rfun.R"); source("rlab_x.R")


# sim params and preparation ----------------------------------------------

p <- 4
N <- 1e3
z <- NULL
for (i in 1:p) z <- cbind(z, arima.sim(list(ar=runif(1,-1,1)),N))
Omega   <- matrix(runif(p^2, 1, 10), ncol = p)
Epsilon <- 1e-4 * matrix(runif(p^2, 1, 10), ncol = p)

  # get simulation data
x <- tvmix(z, Omega, Epsilon)
plot.ts(z)

# perform -----------------------------------------------------------------
dofun <- function(call){
  res <- do.call(call, args = list(x, 12))
  paste(call, MD(res[[1]], Omega))
}

suppressWarnings(
lapply(list("tvsobi2", "SOBI", "tvsobi", "tvsobi011", "tvsobi010", "tvsobi001", "tvsobi000" ), dofun)
)