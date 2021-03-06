# rSim.R contains functions that simulate different type of sources
# rfun.R contains algorithms and utilities of LTV-SOBI
# rlab_x.R is an early work to implement Yeredors' TV-SOBI
# this file runs simulations and save performances to file

remove(list = ls()) # clear environment
library(tidyverse)
source("rfun.R"); source("rsim.R"); source("rlab_x.R");
path <- paste0(getwd(), "/sim/")


# simulation sampling interval --------------------------------------------

do_it_once <- function(x, z, lll = 6, id = "hi", Omega, Epsilon){
  
  for (i in 1:8) {
    flag <- TRUE
    tryCatch({
      if(i == 1) bss_res <- JADE::SOBI(x, k = lll)
      if(i == 2) bss_res <- tvsobi  (x, lag.max = lll, TRUE)
      if(i == 3) bss_res <- tvsobi  (x, lag.max = lll, FALSE)
      if(i == 4) bss_res <- ltvsobi (x, lags = lll, quadratic = TRUE,  fix_symmetry = TRUE,  use_vec = FALSE)
      if(i == 5) bss_res <- ltvsobi (x, lags = lll, quadratic = TRUE,  fix_symmetry = FALSE, use_vec = FALSE)
      if(i == 6) bss_res <- ltvsobi (x, lags = lll, quadratic = FALSE, fix_symmetry = TRUE,  use_vec = FALSE)
      if(i == 7) bss_res <- ltvsobi (x, lags = lll, quadratic = FALSE, fix_symmetry = FALSE, use_vec = FALSE)
      if(i == 8) bss_res <- ltvsobi2(x, lags = lll, use_vec = FALSE)
    }, error = function(e) {
      print(paste("skip to next method due to", e))
      flag <<- FALSE
    })
    
    if(flag) {
      #save_estimator(bss_res, id)
      #save_restored(bss_res, id)
      benchmarks <- SIR_all(bss_res, Omega, Epsilon, z)
      remove(bss_res)
      save_eval(benchmarks, id)
    }
  }  
}


save_eval <- function(benchmarks, id){
  df <- NULL
  for(i in 2:length(benchmarks)){
    df <- data.frame(criteria = attributes(benchmarks)$names[i],
                     value    = benchmarks[[i]]) %>% rbind(df, .)
  }
  df$detail = benchmarks$method
  df$method = word(benchmarks$method)
  df$id     = id
  df$desc   = str_remove(benchmarks$method, " NearestSPD")
  df$N      = benchmarks$N
  df$p      = benchmarks$p
  
  df <- df %>% filter(criteria != "N" & criteria != "p")

  Sys.sleep(runif(1))
  fname <- paste0(path, "benchmarks-", id, ".rds")
  if (file.exists(fname)) saveRDS(rbind(readRDS(fname), df), file = fname)
  else  saveRDS(df, file = fname)
}


# call boosting run -------------------------------------------------------

multido <- function(E, N, sn){
  for(i in 1:50){
    
    Omega   <- matrix(c(2, -9, -4, -6, 5, 6, 0.5, 3, 8), ncol =3)
    Epsilon <- 10^(-E) * matrix(c(-3, -4, 9, 6, 2.5, 2.1, -6, 6, 7), ncol = 3)
    
    zall <- sim_good_sources(10^N, 3, sim_one = "ma", sim_many = "ma")
    xall <- tvmix(zall, Omega, Epsilon)
    
    # loop for freqs
    freq_list <- 2^(0:10)
    for(freq in freq_list){
      for(l in c(3,6,12,1)){
        ids  <- seq(from = 1, to = nrow(xall), by = freq)
        x <- xall[ids,]
        z <- zall[ids,]
        do_it_once(x, z, lll = l, id = paste0("seq", sn, "_ALL_MA_E", E, "N", N, "_Boot_lag", l), Omega, Epsilon)
      }
    }
  }
}
multido(3, 4, 101)
multido(2, 4, 101)

# 
# 
# library(parallel)
# mclapply(
#   as.list(401:432), 
#   function(seq) { 
#     multido(5,5,seq); #multido(5,4,seq);
#     multido(4,5,seq); #multido(4,4,seq);
#   },
#   mc.cores = detectCores()
# )

