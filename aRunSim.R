remove(list = ls()) # clear environment
source("rfun.R"); source("rlab_x.R")


# sim params and preparation ----------------------------------------------

N <- 1e3
z <-cbind(var1 = arima.sim(list(ar=runif(1,-1,1)),N),
          var2 = arima.sim(list(ar=runif(1,-1,1)),N),
          var3 = arima.sim(list(ar=runif(1,-1,1)),N))
Omega <- matrix(runif(9, 1, 10), ncol = 3)
Epsilon <- 1e-4 * matrix(runif(9, 1, 10), ncol = 3)

  # get simulation data
x <- tvmix(z, Omega, Epsilon)


# perform -----------------------------------------------------------------
dofun <- function(call){
  res <- do.call(call, args = list(x, 12))
  paste(call, MD(res[[1]], Omega))
}

lapply(list("SOBI", "tvsobi", "tvsobi0", "tvsobi00", "tvsobi2"), dofun)






