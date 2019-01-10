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
  return(list(X=X, S=z, Omega_true = omega, Epsilon_true = epsilon)) 
}

sim_Generate_2d_Data <- function(N, omega = matrix(c(3, -2, 1, 4), ncol = 2), 
                                 epsilon = matrix(c(-1, -2, 0.5, 1), ncol = 2) *1e-4){
  z1 <- arima.sim(list(ar=c(0.3)),N)
  z2 <- arima.sim(list(ar=c(-0.6)),N)
  z <- cbind(z1,z2)
  X <- matrix(nrow = nrow(z), ncol = ncol(z))
  for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(2) + i * t(epsilon))
  return(list(X=X, S=z, Omega_true = omega, Epsilon_true = epsilon)) 
}

sim_Generate_3d_Data <- function(N, omega = matrix(c(6, -2, 1, -3, -1, 2, 5, 4, 1), ncol = 3),
                                 epsilon = matrix(c(1, 4, 2, -0.8, 3, 0.2, 5, -3, -4), ncol = 3) * 1e-4){
  z1 <- arima.sim(list(ar=c(0.3)),N)
  z2 <- arima.sim(list(ar=c(-0.6)),N)
  z3 <- arima.sim(list(ar=c(0.9)),N)
  z <- apply(cbind(z1,z2,z3), 2, scale)
  X <- matrix(nrow = nrow(z), ncol = ncol(z))
  for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(3) + i * t(epsilon))
  return(list(X=X, S=z, Omega_true = omega, Epsilon_true = epsilon)) 
}


# sim all and measure -----------------------------------------------------

sim_Measure_Performance <- function(X, S, omega, epsilon,
                                    call_method = function(X) tvsobi(X),
                                    sim_id = "manual"){
  N <- nrow(X)
  res <- tryCatch(call_method(X), error = function(e) NA)
  if(!is.na(res)) if(is.null(res$Epsilon)) res$Epsilon <- NA

  md_initial  <- tryCatch( MD(res$W, omega),
                           error = function(e) NA)
  md_middle   <- tryCatch( getMD_t(omega, epsilon, floor(N/2), res$W, res$Epsilon), 
                           error = function(e) NA)
  md_end      <- tryCatch( getMD_t(omega, epsilon, N, res$W, res$Epsilon), 
                           error = function(e) NA)
  md_AVE      <- tryCatch( getMD_ave(omega, epsilon, N, res$W, res$Epsilon), 
                           error = function(e) NA)
  sir_db      <- tryCatch( SIR(res$S, S), 
                           error = function(e) NA) 
  
  
  nearestDist <- tryCatch(ifelse(is.null(res$nearestDist), 0, res$nearestDist), error = function(e) NA)
  
  
  res_W <- tryCatch(ifelse(is.null(res$W), NA, paste(res$W, collapse = ", ")), error = function(e) NA)
  res_Epsilon <- tryCatch(ifelse(is.null(res$Epsilon), NA, paste(res$Epsilon, collapse = ", ")), error = function(e) NA)
  
  data.frame(sim_id = sim_id, N = N, p = ncol(X), method = capture.output(call_method)[1], 
             md_initial = md_initial, md_middle = md_middle, md_end = md_end, md_AVE = md_AVE, sir_db = sir_db,
             nearestFixDist = nearestDist,
             omega_true = paste(omega, collapse = ", "),  epsilon_true = paste(epsilon, collapse = ", "),
             w_result = res_W, epsilon_result = res_Epsilon)
}



sim_All_Mehtods <- function(n_vector = 100 * 2 ^ {0 : 3},
                            sim_funs = list(function(N) sim_Generate_2d_Data(N),
                                             function(N) sim_Generate_3d_Data(N),
                                             function(N) sim_Generate_Yeredor_Data(N)),
                            bss_funs = list(function(X) SOBI(X),
                                             function(X) JADE(X),
                                             # function(X) tvsobi(X, useQuadratic = F, epsilon.method = 1),
                                             # function(X) tvsobi(X, useQuadratic = F, epsilon.method = 2),
                                             # function(X) tvsobi(X, useQuadratic = T, epsilon.method = 1),
                                             # function(X) tvsobi(X, useQuadratic = T, epsilon.method = 2),
                                             function(X) tvsobi(X, useQuadratic = T, epsilon.method = 3)),
                            printProgress = T) {
  
  # init a correct df
  init_sim <- sim_Generate_2d_Data(1e2)
  res_df <- sim_Measure_Performance(init_sim$X, init_sim$S, init_sim$Omega_true, init_sim$Epsilon_true)
  res_df <- res_df[-1,]; remove(init_sim);

  # loop through candidates

  for(N in n_vector){
    for(sim in sim_funs){
      signals <- sim(N)   # result comparable for same N and same sim
      sim_id <- Sys.time() %>% as.double()
      for(bss in bss_funs){
        if(printProgress) cat(N, " |\t ", capture.output(sim)[1], " |\t ", capture.output(bss)[1], "\r")
        res_current <- sim_Measure_Performance(signals$X, signals$S, signals$Omega_true, signals$Epsilon_true,
                                               call_method = bss, as.character(sim_id))
        res_df <- rbind(res_df, res_current)
      }
    }
  }
  
  cat(rep(" ", 80),"\r\n") #clean up print
  return(res_df)
}


quick_check <- function(){
  sim_All_Mehtods(100 * 2^{0:8}) -> aaa
  require(tidyverse)
  aaa %>% 
    ggplot(aes(x = N, y = sir_db)) +
    geom_path()
    facet_wrap(p) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}


sim_bootstrap <- function(n_vector = 100 * 2 ^ {0 : 10}, boot_n = 10, 
                          savefile = "res_PAR_",
                          sqlconn = NA) {
  # init
  res_df <- sim_All_Mehtods(n_vector)
  # boot
  for (i in 2:boot_n) {
    paste(Sys.time(), " runing bootstrap series", i, "/", boot_n) %>% print()
    res_df <- sim_All_Mehtods(n_vector) %>% rbind(res_df)
  }
  
  # save, auto detect the next file to save
  if(!is.na(savefile)){
    seq <- 1 + list.files() %>% 
      str_extract(paste0(savefile, "[0-9][0-9].rds")) %>% 
      str_remove(savefile) %>%
      str_remove(".rds") %>% 
      parse_integer() %>% 
      max(na.rm = T)
    saveRDS(res_df, file = paste0(savefile, seq, ".rds"))
  }
  
  
  
  
}


# CURRENT -----------------------------------------------------------------

hi <- function(){
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
  
}


# parallel ----------------------------------------------------------------

library(parallel)
mclapply(10:50, 
         function(vec) sim_bootstrap(100*2^{0:10}, vec, savefile = "res_PAR_"),
         mc.cores = 15)


df <- readRDS("res_PAR_10.rds")
for(i in 10:50){
  fname <- paste0("res_PAR_", i, ".rds")
  cur <- readRDS(fname)
  df <- rbind(df, cur)
}

file.exists("res_PAR_10.rds")


# lapply(100*2^{0:10}, function(vec) sim_bootstrap(vec, 10, "res_sim_boot.rds"))

# parallel::mclapply()
# map() to boot_n, modify boot_n to vector
