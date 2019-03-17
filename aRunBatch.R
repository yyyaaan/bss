library(bigrquery)
library(parallel)
library(JADE)
remove(list = ls()) # clear environment
source("rfun.R"); source("rlab_x.R")
options(stringsAsFactors = FALSE)
#N <- 1000; p <- 4; id = "TEST"; append = TRUE


# single sim with FIXED source --------------------------------------------

do_one_sim <- function(p = 9, N = 1e3, id = "TEST"){
  
  # get source and mix
  Omega   <- matrix(runif(p^2, 1, 10), ncol = p)
  Epsilon <- 1e-5 * matrix(runif(p^2, 1, 10), ncol = p)
  z <- NULL
  for (i in 1:p){
    rand <- sample(1:3, 1)
    if (rand == 1) var <- arima.sim(list(ar=runif(1,-1,1)),N)
    if (rand == 2) var <- arima.sim(list(ma=runif(1,-1,1)),N)
    if (rand == 3) var <- sin((1:N)/(N/100))
    z   <- cbind(z, var)
  }
  x <- tvmix(z, Omega, Epsilon)
  
  #save source
  #save_signal(z, id, "true_signal")
  #save_signal(z, id, "mixture")
  save_params(Omega, Epsilon, id, "true_param")
  
  # call all methods
  for (i in 1:8) {
    if(i == 1) bss_res <- JADE::SOBI(x, k = 6)
    if(i == 2) bss_res <- tvsobi  (x, lag.max = 6, TRUE)
    if(i == 3) bss_res <- tvsobi  (x, lag.max = 6, FALSE)
    if(i == 4) bss_res <- ltvsobi (x, lags = 6, TRUE, TRUE)
    if(i == 5) bss_res <- ltvsobi (x, lags = 6, TRUE, FALSE)
    if(i == 6) bss_res <- ltvsobi (x, lags = 6, FALSE, TRUE)
    if(i == 7) bss_res <- ltvsobi (x, lags = 6, FALSE, FALSE)
    if(i == 8) bss_res <- ltvsobi2(x, lags = 6)
    
    save_estimator(bss_res, id)
    #save_restored(bss_res, id)
    benchmarks <- SIR_all(bss_res, Omega, Epsilon, z)
    remove(bss_res)
    save_eval(benchmarks, id)
  }  
}



# functions to save to bigquery --------------------------------------------

reset_access_cred(); set_service_token("/home/yanpan/.gcp.json")
bqcon <- dbConnect(bigrquery::bigquery(), project = "yyyaaannn", dataset = "BSS", billing = "yyyaaannn")
# library(odbc); bqcon <- dbConnect(odbc::odbc(), "Study Database")
path <- paste0(getwd(), "/sim/")

save_signal <- function(z, id, type){
  data.frame(id   = id,
             type = type,
             dim  = paste0(nrow(z), "*", ncol(z)),
             t    = 1:nrow(z),
             z    = z) -> df
  colnames(df)[5:ncol(df)] <- paste0("signal_est_", 1:ncol(z))
  Sys.sleep(runif(1))
  saveRDS(df, file = paste0(path, "signals-", id, "-", as.numeric(Sys.time()), ".rds"))
}

save_params <- function(Omega, Epsilon, id, type = "true_param"){
  data.frame(id      = id,
             dim     = paste0(nrow(Omega), "*", ncol(Omega)),
             type    = type,
             omega_vec   = as.vector(Omega),
             epsilon_vec = as.vector(Epsilon)) -> df
  Sys.sleep(runif(1))
  saveRDS(df, file = paste0(path, "params-", id, "-", as.numeric(Sys.time()), ".rds"))
  # tryCatch(dbWriteTable(bqcon, "params", df, append = TRUE), 
  #          error = function(e) {saveRDS(df, file = paste0(path, "params-", id, "-", as.numeric(Sys.time()), ".rds")); print(e)})
}

save_eval <- function(benchmarks, id){
  df <- NULL
  for(i in 2:length(benchmarks)){
    df <- data.frame(criteria = attributes(benchmarks)$names[i],
                     value    = benchmarks[[i]]) %>% rbind(df, .)
  }
  df$method = benchmarks$method
  df$id     = id
  Sys.sleep(runif(1))
  saveRDS(df, file = paste0(path, "benchmarks-", id, "-", as.numeric(Sys.time()), ".rds"))
}
  

save_estimator <- function(bss_res, id){
  if (class(bss_res)=="bss") {
    Omega    <- solve(bss_res$W)
    Epsilon <- 0
  }
  if (class(bss_res)== "tvbss"){
    Omega   <- bss_res$Omega_hat
    Epsilon <- bss_res$Epsilon_hat
  }
  
  save_params(Omega, Epsilon, id, type = bss_res$method)

}
save_restored <- function(bss_res, id){
  save_signal(z = bss_res$S, id = id, type = bss_res$method)
}




# function for Parallel Computing -----------------------------------------

do_full_set <- function(rid = "XXXX"){
  Ns <- 50 * 2 ^ {1:11}
  ps  <- 3:9
  for(N in Ns) for(p in ps) {
    tryCatch(do_one_sim(p, N, paste0( as.numeric(Sys.Date()), "ID", rid, "D", p, "N", N)),
             error = function(e) print(paste("error occured and ignored:", e)))
  }
  print(paste("BATCH", rid, "Completed"))
}

for(i in 1:30){
  
}

for (i in 1:100) do_full_set(9000 + i)
# mclapply(9001:9399, function(rid) do_full_set(rid), mc.cores = detectCores())


# file processing ---------------------------------------------------------

combined_files <- function(indicator = "benchmarks"){
  all_files <- list.files("./sim")
  all_files <- all_files[grep(".rds", all_files)]
  all_files <- all_files[grep(indicator, all_files)]
  df <- NA
  for(f in all_files){
    read_df <- readRDS(paste0("./sim/", f))  
    df <- rbind(df, read_df)
  }
  saveRDS(df, file = paste0(indicator, ".rds"))
}
