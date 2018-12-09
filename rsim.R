source("rlab.R")
library(magrittr)
library(ggplot2)

# new simulation ----------------------------------------------------------

sim_Generate_2d_Data <- function(N, omega = matrix(c(3, -2, 1, 4), ncol = 2), 
                                 epsilon = matrix(c(-1, -2, 0.5, 1), ncol = 2) *1e-4){
  # yeredor's simulation
  z <- rnorm(N+4)
  z1 <- 1 + 2*z[1:N + 3] - 0.5*z[1:N + 2] - z[1:N + 1] + z[1:N]
  z <- rnorm(N+3)
  z2 <- 1 - z[1:N + 2] + 3*z[1:N + 1] + 2*z[1:N]
  z <- cbind(z1,z2)
  X <- matrix(nrow = nrow(z), ncol = ncol(z))
  for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(2) + i * t(epsilon))
  return(X) 
}

sim_Generate_3d_Data <- function(N, omega = matrix(c(6, -2, 1, -3, -1, 2, 5, 4, 1), ncol = 3),
                                 epsilon = matrix(c(1, 4, 2, -0.8, 3, 0.2, 5, -3, -4), ncol = 3) * 1e-4){
  z1 <- arima.sim(list(ar=c(0.3)),N)
  z2 <- arima.sim(list(ar=c(-0.6)),N)
  z3 <- arima.sim(list(ar=c(0.9)),N)
  z <- apply(cbind(z1,z2,z3), 2, scale)
  X <- matrix(nrow = nrow(z), ncol = ncol(z))
  for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(3) + i * t(epsilon))
  return(X)
}

sim_Measure_Performance <- function(X, omega, epsilon,
                                    msg = "TVSOBI, Quadratic TRUE, Epsilon Mehtod 1",
                                    call_method = function(X) tvsobi(X)){
  N <- nrow(X)
  start_time <- Sys.time()
  res <- tryCatch(call_method(X), error = function(e) NA)
  md_initial <- tryCatch(MD(res$W, omega), error = function(e) NA)
  md_middle <- tryCatch(getMD_t(res, omega, epsilon, floor(N/2)), error = function(e) NA)
  md_end <- tryCatch(getMD_t(res, omega, epsilon, N), error = function(e) NA)
  nearestDist <- tryCatch(ifelse(is.null(res$nearestDist), 0, res$nearestDist), error = function(e) NA)
  time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  data.frame(N = N, p = ncol(X), method = msg, 
             md_initial = md_initial, md_middle = md_middle, md_end = md_end,
             time = time, nearestFixDist = nearestDist,
             omega = paste(omega, collapse = ", "),
             epsilon = paste(epsilon, collapse = ", "))
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
        for(i_opt in c(1,2,3)){
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

# current lab -------------------------------------------------------------

# sim_single <- sim_All(n_vector = 100 * 2 ^ {0 : 9})
# myColor <- RColorBrewer::brewer.pal(3, "Set2")
# ggplot(sim_single) +
#   geom_line(aes(N, md_initial), color = myColor[1]) +
#   geom_line(aes(N, md_middle), color = myColor[2]) +
#   geom_line(aes(N, md_end), color = myColor[3]) +
#   facet_wrap(~method)

temp <- readRDS("res_sim_boot.rds"); temp <- temp[-(1:nrow(temp)), ]; saveRDS(temp, "res_sim_boot.rds")
sim_bootstrap(100*2^{0:10}, 99, "res_sim_boot.rds")

# call_method = function(x) tvsobi(x); msg = "test";
# X = sim_Generate_2d_Data(1e2);  omega = matrix(c(3, -2, 1, 4), ncol = 2); epsilon = matrix(c(-1, -2, 0.5, 1), ncol = 2) *1e-4;
# X = sim_Generate_3d_Data(1e2); omega = matrix(c(6, -2, 1, -3, -1, 2, 5, 4, 1), ncol = 3);  epsilon = matrix(c(1, 4, 2, -0.8, 3, 0.2, 5, -3, -4), ncol = 3) * 1e-4
