source("rlab.R")
library(magrittr); library(JADE)
options(stringsAsFactors = F)

# orignal parameters
# omega_2d <- matrix(c(3, -2, 1, 4), ncol = 2)
# epsilon_2d <- matrix(c(-1, -2, 0.5, 1), ncol = 2) *1e-4
# omega_3d <- matrix(c(6, -2, 1, -3, -1, 2, 5, 4, 1), ncol = 3)
# epsilon_3d <- matrix(c(1, 4, 2, -0.8, 3, 0.2, 5, -3, -4), ncol = 3) * 1e-4


# new simulation ----------------------------------------------------------

sim_Generate_Yeredor_Data <- function(N, omega = matrix(runif(4,-10,10), ncol=2), 
                                      epsilon = matrix(runif(4,-10,10), ncol=2)*1e-4,
                                      info = "rnd ar-rnd"){
  # yeredor's simulation
  z <- rnorm(N+4)
  z1 <- 1 + 2*z[1:N + 3] - 0.5*z[1:N + 2] - z[1:N + 1] + z[1:N]
  z <- rnorm(N+3)
  z2 <- 1 - z[1:N + 2] + 3*z[1:N + 1] + 2*z[1:N]
  z <- cbind(z1,z2)
  return(list(X=make_tvmix(z, omega, epsilon), S=z, Omega_true = omega, Epsilon_true = epsilon)) 
}

sim_Generate_2d_Data <- function(N, omega = matrix(runif(4,-10,10), ncol=2), 
                                 epsilon = matrix(runif(4,-10,10), ncol=2)*1e-4){
  z1 <- arima.sim(list(ar=runif(1,-1,1)),N)
  z2 <- arima.sim(list(ar=runif(1,-1,1)),N)
  z <- cbind(z1,z2)
  X <- matrix(nrow = nrow(z), ncol = ncol(z))
  for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(2) + i * t(epsilon))
  return(list(X=X, S=z, Omega_true = omega, Epsilon_true = epsilon)) 
}

sim_Generate_3d_Data <- function(N, omega = matrix(runif(9,-10,10), ncol=3), 
                                 epsilon = matrix(runif(9,-10,10), ncol=3)*1e-4){
  z1 <- arima.sim(list(ar=runif(1,-1,1)),N)
  z2 <- arima.sim(list(ar=runif(1,-1,1)),N)
  z3 <- arima.sim(list(ar=runif(1,-1,1)),N)
  z <- apply(cbind(z1,z2,z3), 2, scale)
  X <- matrix(nrow = nrow(z), ncol = ncol(z))
  for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(3) + i * t(epsilon))
  return(list(X=X, S=z, Omega_true = omega, Epsilon_true = epsilon)) 
}

sim_Generate_Random_Data <- function(N, dim = 4, level = 1e-4){
  p <- dim
  omega   <- matrix(runif(p^2,-10,10), ncol=p)
  epsilon <- matrix(runif(p^2,-10,10), ncol=p)*level
  z <- matrix(nrow = N, ncol = p)
  for (i in 1:p) z[,i] <- arima.sim(list(ar = runif(1, -1, 1)), N)
  z <- apply(z, 2, scale)
  X <- matrix(nrow = nrow(z), ncol = ncol(z))
  for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(p) + i * t(epsilon))
  return(list(X=X, S=z, Omega_true = omega, Epsilon_true = epsilon)) 
}


# sim all and measure -----------------------------------------------------

sim_Measure_Performance <- function(X, S, omega, epsilon,
                                    call_method = function(X) tvsobi(X),
                                    sim_id = "manual", sim_fun = "unknown"){
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
  
  nearestDist <- tryCatch(ifelse(is.null(res$nearestDist), 0, res$nearestDist), 
                          error = function(e) NA)
  
  res_W       <- tryCatch(ifelse(is.null(res$W), NA, paste(res$W, collapse = ", ")), 
                          error = function(e) NA)
  res_Epsilon <- tryCatch(ifelse(is.null(res$Epsilon), NA, paste(res$Epsilon, collapse = ", ")), 
                          error = function(e) NA)
  
  data.frame(sim_id = sim_id, N = N, p = ncol(X), 
             sim_fun = sim_fun, bss_fun = capture.output(call_method)[1], 
             md_initial = md_initial, md_middle = md_middle, md_end = md_end, md_AVE = md_AVE, sir_db = sir_db,
             nearestFixDist = nearestDist,
             omega_true = paste(omega, collapse = ", "),  epsilon_true = paste(epsilon, collapse = ", "),
             w_result = res_W, epsilon_result = res_Epsilon)
}



sim_All_Mehtods <- function(n_vector = 100 * 2 ^ {0 : 3},
                            sim_funs = list(function(N) sim_Generate_Random_Data(N, dim=4, level=1e-5),
                                            function(N) sim_Generate_Random_Data(N, dim=5, level=1e-5),
                                            function(N) sim_Generate_Random_Data(N, dim=6, level=1e-5),
                                            function(N) sim_Generate_Random_Data(N, dim=7, level=1e-5),
                                            function(N) sim_Generate_Random_Data(N, dim=8, level=1e-5),
                                            function(N) sim_Generate_Random_Data(N, dim=9, level=1e-5)),
                            bss_funs = list(function(X) SOBI(X),
                                            function(X) JADE(X),
                                            # function(X) tvsobi(X, useQuadratic = F, epsilon.method = 1),
                                            # function(X) tvsobi(X, useQuadratic = F, epsilon.method = 2),
                                            # function(X) tvsobi(X, useQuadratic = T, epsilon.method = 1),
                                            # function(X) tvsobi(X, useQuadratic = T, epsilon.method = 2),
                                            # function(X) tvsobi(X, useQuadratic = T, epsilon.method = 4),
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
                                               call_method = bss, 
                                               as.character(sim_id), capture.output(sim)[1])
        res_df <- rbind(res_df, res_current)
      }
    }
  }
  
  if(printProgress) cat(rep(" ", 200),"\r\n") #clean up print
  return(res_df)
}


sim_bootstrap <- function(n_vector = 100 * 2 ^ {0 : 10}, boot_n = 10, 
                          savefile = "res_PAR_",
                          sqlconn = NA, sqltable = NA) {
  # init data frame
  if(length(n_vector) * boot_n > 999) warning("large bootstrap, consider SQL or lapply to avoid crash")
  res_df <- sim_All_Mehtods(n_vector)
  
  #init sql table
  if(!is.na(sqlconn)) sqlconn %>% dbWriteTable(sqltable, res_df, append = dbExistsTable(., sqltable)) 
  
  
  # boot and save to sql
  for (i in 2:boot_n) {
    
    res_cur <- sim_All_Mehtods(n_vector, printProgress = F)
    
    if(!is.na(savefile)){
      # without SQL, all result will be kept in system RAM
      res_df <- rbind(res_df, res_cur)
    } 
    
    if(!is.na(sqlconn)){
      # with SQL, append rows efficiently
      sqlconn %>% dbWriteTable(sqltable, res_cur, append = T) #init sql table 
      print(paste("append to SQL ok", i, "of", boot_n))
    }
  }
  
  # save, auto detect the next file to save
  if(!is.na(savefile)){
    require(tidyverse)
    seq <- 1 + list.files() %>% 
      str_extract(paste0(savefile, "[0-9][0-9].rds")) %>% 
      str_remove(savefile) %>%
      str_remove(".rds") %>% 
      parse_integer() %>% 
      max(na.rm = T)
    if(seq < 0) seq <- 10
    saveRDS(res_df, file = paste0(savefile, seq, ".rds"))
    paste(Sys.time(), " current bootstrap of", boot_n, "compelted successfully | file saved:", seq) %>% print()
  }
  
  # save to sql, append, increase performance
}


# quick check -------------------------------------------------------------

quick_check <- function(){
  require(tidyverse)
  
  aaa <- sim_All_Mehtods(
    n_vector = 100 * 2 ^ {0 : 5},
    sim_funs = list(function(N) sim_Generate_2d_Data(N, matrix(runif(4,-10,10), ncol=2), matrix(runif(4,-10,10), ncol=2)*1e-3),
                    function(N) sim_Generate_3d_Data(N, matrix(runif(9,-10,10), ncol=3), matrix(runif(9,-10,10), ncol=3)*1e-3),
                    function(N) sim_Generate_Yeredor_Data(N, matrix(runif(4,-10,10), ncol=2), matrix(runif(4,-10,10), ncol=2)*1e-3)),
    bss_funs = list(function(X) SOBI(X),
                    function(X) JADE(X),
                    function(X) tvsobi(X, useQuadratic = F, epsilon.method = 1),
                    function(X) tvsobi(X, useQuadratic = F, epsilon.method = 2),
                    function(X) tvsobi(X, useQuadratic = T, epsilon.method = 1),
                    function(X) tvsobi(X, useQuadratic = T, epsilon.method = 2),
                    function(X) tvsobi(X, useQuadratic = T, epsilon.method = 4),
                    function(X) tvsobi(X, useQuadratic = T, epsilon.method = 3)),
    printProgress = T)
  
  
  aaa %>% 
    gather("key", "value", md_AVE, md_initial, sir_db) %>%
    ggplot(aes(x = N, y = value, color = bss_fun)) +
    geom_path() +
    facet_grid(key~sim_fun, scales="free_y") +
    scale_x_log10(limits = c(min(aaa$N), max(aaa$N))) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "bottom")
}


# parallel and save file --------------------------------------------------


parallel_save <- function(){
  library(parallel)
  mclapply(10:50,
           function(vec) sim_bootstrap(100*2^{0:10}, vec, savefile = "res_PAR_"),
           mc.cores = detectCores() - 1)
}

combine_PARs <- function(){
  res <- readRDS("res_PAR_10.rds")
  for (ii in 11:50) res <- rbind(res, readRDS(paste0("res_PAR_", ii, ".rds")))
  saveRDS(res, file = "res_boot_sir_1000.rds")
}



# SQL parallel ------------------------------------------------------------

require(odbc); require(parallel)

sqlconn <- dbConnect(odbc::odbc(), "Study Database")
mclapply(rep(8,120),
         function(x) sim_bootstrap(100*2^{0:10}, x, savefile = NA,
                                   sqlconn = sqlconn, sqltable = "boot_rnd-5_rndAR_190115"),
         mc.cores = detectCores() - 1)