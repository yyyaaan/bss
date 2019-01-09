source("rlab.R")
library(magrittr); library(JADE)
options(stringsAsFactors = F)

# new simulation ----------------------------------------------------------

sim_Generate_Yeredor_Data <- function(N, omega = matrix(c(3, -2, 1, 4), ncol = 2), 
                                 epsilon = matrix(c(-1, -2, 0.5, 1), ncol = 2) *1e-4){
  # yeredor's simulation
  z <- rnorm(N+4)
  z1 <- 1 + 2*z[1:N + 3] - 0.5*z[1:N + 2] - z[1:N + 1] + z[1:N]
  z <- rnorm(N+3)
  z2 <- 1 - z[1:N + 2] + 3*z[1:N + 1] + 2*z[1:N]
  z <- cbind(z1,z2)
  X <- matrix(nrow = nrow(z), ncol = ncol(z))
  for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(2) + i * t(epsilon))
  return(list(X=X, S=z)) 
}

sim_Generate_2d_Data <- function(N, omega = matrix(c(3, -2, 1, 4), ncol = 2), 
                                 epsilon = matrix(c(-1, -2, 0.5, 1), ncol = 2) *1e-4){
  z1 <- arima.sim(list(ar=c(0.3)),N)
  z2 <- arima.sim(list(ar=c(-0.6)),N)
  z <- cbind(z1,z2)
  X <- matrix(nrow = nrow(z), ncol = ncol(z))
  for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(2) + i * t(epsilon))
  return(list(X=X, S=z)) 
}

sim_Generate_3d_Data <- function(N, omega = matrix(c(6, -2, 1, -3, -1, 2, 5, 4, 1), ncol = 3),
                                 epsilon = matrix(c(1, 4, 2, -0.8, 3, 0.2, 5, -3, -4), ncol = 3) * 1e-4){
  z1 <- arima.sim(list(ar=c(0.3)),N)
  z2 <- arima.sim(list(ar=c(-0.6)),N)
  z3 <- arima.sim(list(ar=c(0.9)),N)
  z <- apply(cbind(z1,z2,z3), 2, scale)
  X <- matrix(nrow = nrow(z), ncol = ncol(z))
  for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(3) + i * t(epsilon))
  return(list(X=X, S=z)) 
}


# sim all and measure -----------------------------------------------------

sim_Measure_Performance <- function(X, omega, epsilon,
                                    msg = "TVSOBI, Quadratic TRUE, Epsilon Mehtod 1",
                                    call_method = function(X) tvsobi(X)){
  N <- nrow(X)
  start_time <- Sys.time()
  res <- tryCatch(call_method(X), error = function(e) NA)
  time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  md_initial <- tryCatch(MD(res$W, omega), error = function(e) NA)
  md_middle <- tryCatch(getMD_t(res, omega, epsilon, floor(N/2)), error = function(e) NA)
  md_end <- tryCatch(getMD_t(res, omega, epsilon, N), error = function(e) NA)
  md_AVE <- tryCatch(getMD_ave(res, omega, epsilon, N), error = function(e) NA)
  nearestDist <- tryCatch(ifelse(is.null(res$nearestDist), 0, res$nearestDist), error = function(e) NA)
  res_W <- tryCatch(ifelse(is.null(res$W), NA, paste(res$W, collapse = ", ")), error = function(e) NA)
  res_Epsilon <- tryCatch(ifelse(is.null(res$Epsilon), NA, paste(res$Epsilon, collapse = ", ")), error = function(e) NA)
  
  data.frame(N = N, p = ncol(X), method = msg, 
             md_initial = md_initial, md_middle = md_middle,md_end = md_end, md_AVE = md_AVE,
             time = time, nearestFixDist = nearestDist,
             omega = paste(omega, collapse = ", "),  epsilon = paste(epsilon, collapse = ", "),
             res_w = res_W, res_Epsilon = res_Epsilon)
}

sim_All <- function(n_vector = 100 * 2 ^ {0 : 3} ) {
  # init
  sim_df <- sim_Measure_Performance(X = sim_Generate_2d_Data(1e2), 
                                    omega = matrix(c(3, -2, 1, 4), ncol = 2), 
                                    epsilon = matrix(c(-1, -2, 0.5, 1), ncol = 2) *1e-4)
  sim_df <- sim_df[-1,]

  # params i_x = 1,2
  epsilon <- list(matrix(c(1, 4, 2, -0.8, 3, 0.2, 5, -3, -4), ncol = 3) * 1e-4,
                  matrix(c(-1, -2, 0.5, 1), ncol = 2) *1e-4)
  omega <- list(matrix(c(6, -2, 1, -3, -1, 2, 5, 4, 1), ncol = 3),
                matrix(c(3, -2, 1, 4), ncol = 2))


  
  # loop
  for (n in n_vector) {
    paste(Sys.time(), "runing with n =", n) %>% print()
    for (i_x in 1:2) {
      epsilon_cur <- epsilon[[i_x]]
      omega_cur <- omega[[i_x]]
      if(dim(epsilon_cur)[1] == 2) X <- sim_Generate_2d_Data(n, omega_cur, epsilon_cur)
      if(dim(epsilon_cur)[1] == 3) X <- sim_Generate_3d_Data(n, omega_cur, epsilon_cur)

      # trying different methods
      sim_Measure_Performance(X, omega_cur, epsilon_cur, 
                              call_method = function(X) SOBI(X),
                              msg = "SOBI") %>%
        rbind(sim_df) -> sim_df
      for (i_quad in c(TRUE, FALSE)) {
  ### no longer simulate for method 3
        for(i_opt in c(1,2)){
          sim_Measure_Performance(X, omega_cur, epsilon_cur, 
                                  call_method = function(X) tvsobi(X, useQuadratic = i_quad, epsilon.method = i_opt),
                                  msg = paste("TVSOBI, Quadratic", i_quad, ", Epsilon Mehtod", i_opt)) %>%
            rbind(sim_df) -> sim_df
        }
      }
    }
  }
  
  sim_df
}

sim_bootstrap <- function(n_vector = 100 * 2 ^ {0 : 3}, boot_n = 10, filename) {
  # init
  sim_df <- readRDS(filename)
  for (i in 1:boot_n) {
    cat("runing bootstrap series", i, "/", boot_n, "\n")
    sim_All(n_vector) %>% rbind(sim_df) -> sim_df
    saveRDS(sim_df, filename)
  }
}



# CURRENT -----------------------------------------------------------------

mat0 <- matrix(rep(0,4), ncol = 2)

omega2d = matrix(c(3, -2, 1, 4), ncol = 2)
epsilon2d = matrix(c(-1, -2, 0.5, 1), ncol = 2) *1e-4

aaa <- sim_Generate_2d_Data(1e3, omega2d, epsilon2d)
res_sobi <- SOBI(aaa$X)
res_tvsobi <- tvsobi(aaa$X, epsilon.method = 3)
SIR(aaa$S, res_sobi$S); MD(omega, res_sobi$W)

SIR(aaa$S, res_tvsobi$S); MD(omega, res_tvsobi$W); 
getMD_ave(omega2d, epsilon2d, ncol(aaa$X), res_tvsobi$W, res_tvsobi$Epsilon)


S <- cbind(rt(1000, 4), rnorm(1000), runif(1000))
A <- matrix(rnorm(9), ncol = 3)
X <- S %*% t(A)
SIR(S, JADE(X)$S)




# parallel ----------------------------------------------------------------


# temp <- sim_All(); temp <- temp[-(1:nrow(temp)), ]; saveRDS(temp, "res_sim_boot.rds")

lapply(100*2^{0:10}, function(vec) sim_bootstrap(vec, 10, "res_sim_boot.rds"))

# parallel::mclapply()
# map() to boot_n, modify boot_n to vector

sim_bootstrap(100*2^{0:10}, 1000, "res_sim_boot.rds")



library(parallel)

