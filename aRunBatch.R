library(bigrquery)
library(parallel)
library(JADE)
remove(list = ls()) # clear environment
source("rfun.R"); source("rlab_x.R")
source("rsim.R")
options(stringsAsFactors = FALSE)
#N <- 10000; p <- 4; id = "TEST"; append = TRUE


# single sim with FIXED source --------------------------------------------

do_one_sim <- function(p = 9, N = 1e3, id = "TEST"){
  
  # get source and mix
  Omega   <- matrix(runif(p^2, 1, 10), ncol = p)
  Epsilon <- 1e-4 * matrix(runif(p^2, 1, 10), ncol = p)
  # z <- NULL
  # for (i in 1:p){
  #   rand <- sample(1:3, 1)
  #   if (rand == 1) var <- arima.sim(list(ar=runif(1,-1,1)),N)
  #   if (rand == 2) var <- arima.sim(list(ma=runif(1,-1,1)),N)
  #   if (rand == 3) var <- sin((1:N)/(N/100))
  #   z   <- cbind(z, var)
  # }
  z <- sim_good_sources(N, p)
  x <- tvmix(z, Omega, Epsilon)
  x <- scale(x, scale = FALSE)
  #save source
  #save_signal(z, id, "true_signal")
  #save_signal(z, id, "mixture")
  save_params(Omega, Epsilon, id, "true_param")
  
  # call all methods
  for (i in 1:8) {
    flag <- TRUE
    tryCatch({
      if(i == 1) bss_res <- JADE::SOBI(x, k = 6)
      if(i == 2) bss_res <- tvsobi  (x, lag.max = 6, TRUE)
      if(i == 3) bss_res <- tvsobi  (x, lag.max = 6, FALSE)
      if(i == 4) bss_res <- ltvsobi (x, lags = 6, quadratic = TRUE,  fix_symmetry = TRUE,  use_vec = FALSE)
      if(i == 5) bss_res <- ltvsobi (x, lags = 6, quadratic = TRUE,  fix_symmetry = FALSE, use_vec = FALSE)
      if(i == 6) bss_res <- ltvsobi (x, lags = 6, quadratic = FALSE, fix_symmetry = TRUE,  use_vec = FALSE)
      if(i == 7) bss_res <- ltvsobi (x, lags = 6, quadratic = FALSE, fix_symmetry = FALSE, use_vec = FALSE)
      if(i == 8) bss_res <- ltvsobi2(x, lags = 6, use_vec = FALSE)
    }, error = function(e) {
      print(paste("skip to next due to", e))
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
  Ns <- 50 * 2 ^ {11:1}
  ps  <- 9:3
  for(p in ps) for(N in Ns) {
    tryCatch(do_one_sim(p, N, paste0( as.numeric(Sys.Date()), "ID", rid, "D", p, "N", N)),
             error = function(e) print(paste("error occured and ignored:", e)))
  }
  print(paste("BATCH", rid, "Completed"))
}


for (i in 1:99) do_full_set(paste0("NEW", 6000 + i))
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

# benchmarking plotting ---------------------------------------------------
plotting <- function() {
  library(tidyverse)
  combined_files()
  benchmarks <- readRDS("benchmarks.rds")
  benchmarks$m = str_remove(benchmarks$method, " NearestSPD")
  benchmarks$simple_method = word(benchmarks$method)
  benchmarks$N = str_extract(benchmarks$id, "N[:digit:]{3,8}") %>% str_remove("N") %>% parse_integer()
  benchmarks$p = str_extract(benchmarks$id, "D[:digit:]N") %>% str_remove("D") %>% str_remove("N") %>% parse_integer()
  benchmarks$logN = log(benchmarks$N)
  
  benchmarks %>%
    filter(!is.na(value)) %>%
    filter(p <= 6) %>%
    group_by(criteria, simple_method, N, p, logN) %>%
    summarise_at("value", mean) %>%
    ggplot(aes(x = logN, y = value, color = simple_method)) +
    geom_line() +
    facet_grid(criteria ~ p, scales = "free_y")
  
}
print("Compeleted!")