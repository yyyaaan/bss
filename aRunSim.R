remove(list = ls()) # clear environment
source("rfun.R"); source("rlab_x.R")


# sim params and preparation ----------------------------------------------

p <- 4
N <- 1e3
z <- NULL
Omega   <- matrix(runif(p^2, 1, 10), ncol = p)
Matrix_0 <- matrix(rep(0, p^2), ncol = p)
Epsilon <- 1e-5 * matrix(runif(p^2, 1, 10), ncol = p)


# sim signal type ---------------------------------------------------------

for (i in 1:p){
  rand <- sample(1:3, 1)
  if (rand == 1) var <- arima.sim(list(ar=runif(1,-1,1)),N)
  if (rand == 2) var <- arima.sim(list(ma=runif(1,-1,1)),N)
  if (rand == 3) var <- sin((1:N)/(N/100))
  z   <- cbind(z, var)
}

  # get simulation data
x <- tvmix(z, Omega, Epsilon)
# plot.ts(z); plot.ts(x)


# playground --------------------------------------------------------------

res2 <- JADE::SOBI(x)
res1 <- tvsobi011(x)
res0 <- tvsobi(x)
SIR_all(res0, Omega, Epsilon, z)
SIR_all(res1, Omega, Epsilon, z)
SIR_all(res2, Omega, Epsilon, z)
