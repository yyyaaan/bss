# simulate and evaluate ---------------------------------------------------

  ## params-simulation
N <- 1e3
z <-cbind(var1 = arima.sim(list(ar=runif(1,-1,1)),N),
          var2 = arima.sim(list(ar=runif(1,-1,1)),N),
          var3 = arima.sim(list(ar=runif(1,-1,1)),N))
Omega <- matrix(runif(9, 1, 10), ncol = 3)
Epsilon <- 1e-4 * matrix(runif(9, 1, 10), ncol = 3)
#Epsilon <- (Epsilon + t(Epsilon)) / 2

  ## get simulation data
x <- tvmix(z, Omega, Epsilon)

  ## params-evaluation
lags <- 6
  

tvsobi_sym(x, lags) %>% MD(Omega) %>% print
tvsobi(x, lag.max = lags) %>% use_series(W) %>% MD(Omega) %>% print

x %>% cov_sep %>% use_series(beta_3) %>% 
  approxJD() %>% use_series(W) %>% 
  MD(Omega %*% Epsilon) %>% print

x %>% cov_sep_vec %>% use_series(beta_3) %>% 
  approxJD() %>% use_series(W) %>% 
  MD(Omega %*% Epsilon) %>% print

ss <- tvsobi(x, lag.max = 12); MD(ss$W, Omega); MD(ss$Epsilon, solve(Epsilon))
ss <- tvsobi2(x,lags = 12); MD(ss$Omega_hat %>% solve, Omega);  MD(ss$Epsilon_hat %>% solve, Epsilon)
ss <- tvsobi3(x,lags = 12); MD(ss$Omega_hat %>% solve, Omega);  MD(ss$Epsilon_hat %>% solve, Epsilon)

SIR(z, tvunmix(x, ss$Omega_hat, ss$Epsilon_hat))



# batch sim ---------------------------------------------------------------

require(parallel)
parallel_save <- function(){
  mclapply(10:50,
           function(vec) sim_bootstrap(100*2^{0:10}, vec, savefile = "res_PAR_"),
           mc.cores = detectCores() - 1)
}
