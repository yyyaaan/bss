library(parallel)
library(JADE)
library(tidyverse)
remove(list = ls()) # clear environment
source("rfun.R"); source("rlab_x.R")
source("rsim.R")
# source("bqhelper.R")
options(stringsAsFactors = FALSE)
#N <- 10000; p <- 4; id = "TEST"; append = TRUE


# single sim with FIXED source --------------------------------------------

do_one_sim <- function(N = 1e3, p = 9, id = "TEST"){
  
  # get source and mix
  Omega   <- matrix(runif(p^2, 1, 10), ncol = p)
  Epsilon <- 1e-4 * matrix(runif(p^2, 1, 10), ncol = p)

  z <- sim_good_sources(N, p, sim_one = c("ar", "ma"), sim_many = c("ar", "ma"))
  x <- tvmix(z, Omega, Epsilon)
  x <- scale(x, scale = FALSE)
  
  #save_signal(z, id, "true_signal")
  #save_signal(z, id, "mixture")
  save_params(Omega, Epsilon, id, "true_param")
  
  # call all methods
  for (i in 1:8) {
    flag <- TRUE
    tryCatch({
      if(i == 1) bss_res <- JADE::SOBI(x, k = 12)
      if(i == 2) bss_res <- tvsobi  (x, lag.max = 1, TRUE)
      if(i == 3) bss_res <- tvsobi  (x, lag.max = 1, FALSE)
      if(i == 4) bss_res <- ltvsobi (x, lags = 1, quadratic = TRUE,  fix_symmetry = TRUE,  use_vec = FALSE)
      if(i == 5) bss_res <- ltvsobi (x, lags = 1, quadratic = TRUE,  fix_symmetry = FALSE, use_vec = FALSE)
      if(i == 6) bss_res <- ltvsobi (x, lags = 1, quadratic = FALSE, fix_symmetry = TRUE,  use_vec = FALSE)
      if(i == 7) bss_res <- ltvsobi (x, lags = 1, quadratic = FALSE, fix_symmetry = FALSE, use_vec = FALSE)
      if(i == 8) bss_res <- ltvsobi2(x, lags = 1, use_vec = FALSE)
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



# functions to save to bigquery --------------------------------------------

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
  #save2bq(df, "signals", "BSS", hold = 10)
}

save_params <- function(Omega, Epsilon, id, type = "true_param"){
  data.frame(id      = id,
             dim     = paste0(nrow(Omega), "*", ncol(Omega)),
             type    = type,
             omega_vec   = as.vector(Omega),
             epsilon_vec = as.vector(Epsilon)) -> df
  Sys.sleep(runif(1))
  saveRDS(df, file = paste0(path, "params-", id, "-", as.numeric(Sys.time()), ".rds"))
  #save2bq(df, "params", "BSS", hold = 3)
  
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
  df$N = str_extract(id, "N[:digit:]{3,8}") %>% str_remove("N") %>% parse_integer()
  df$p = str_extract(id, "D[:digit:]N") %>% str_remove("D") %>% str_remove("N") %>% parse_integer()
  
  Sys.sleep(runif(1))
  saveRDS(df, file = paste0(path, "benchmarks-", id, "-", as.numeric(Sys.time()), ".rds"))
  #save2bq(df, "benchmarks", "BSS", hold = 50)
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





# submit batch tasks ------------------------------------------------------

do_full_set <- function(rid = "XXXX"){
  cat("set started", rid, "\n\n\n")
  Ns <- 50 * 2 ^ {11:1}
  ps  <- 9:3
  for(p in ps) for(N in Ns) {
    tryCatch(do_one_sim(N, p, paste0( as.numeric(Sys.Date()), "ID", rid, "D", p, "N", N)),
             error = function(e) print(paste("error occured and ignored:", e)))
  }
  print(paste("BATCH", rid, "Completed"))
}

mclapply(999001:999128, function(rid) do_full_set(paste0("ARMA", rid)), mc.cores = detectCores()-1)








# file processing ---------------------------------------------------------

combine_files <- function(indicator = "benchmarks", process = TRUE){
  library(tidyverse)
  all_files <- list.files("./sim")
  all_files <- all_files[grep(".rds", all_files)]
  all_files <- all_files[grep(indicator, all_files)]
  df <- NA
  for(f in all_files){
    read_df <- readRDS(paste0("./sim/", f))  
    df <- rbind(df, read_df)
  }
  
  # processing
  if(process){
    df$m = str_remove(df$method, " NearestSPD")
    df$simple_method = word(df$method)
    df$N = str_extract(df$id, "N[:digit:]{3,8}") %>% str_remove("N") %>% parse_integer()
    df$p = str_extract(df$id, "D[:digit:]N") %>% str_remove("D") %>% str_remove("N") %>% parse_integer()
    df$logN = log(df$N)
  }
  
  saveRDS(df, file = paste0(indicator, ".rds"))
  df
}

# benchmarking plotting ---------------------------------------------------
plotting <- function() {
  
  combine_files() %>%
    filter(!is.na(value)) %>%
    filter(p <= 6) %>%
    group_by(criteria, simple_method, N, p, logN) %>%
    summarise_at("value", mean) %>%
    ggplot(aes(x = logN, y = value, color = simple_method)) +
    geom_line() +
    facet_grid(criteria ~ p, scales = "free_y")
  
}

for (i in 1:3) do_full_set(paste0("BQBQ", 3000 + i))
# mclapply(9001:9399, function(rid) do_full_set(rid), mc.cores = detectCores())

# progress?
list.files("./sim") %>% substr(1,10)%>% str_extract("NEW[:digit:]{4}") %>% str_remove("NEW") %>% parse_integer() %>% max

print("Compeleted!")








simdata <- readRDS("marks.rds")
temp <- readRDS("xxx.rds")
temp$sim <- "ARMA"
all <- rbind(simdata, temp)
all$detail <- str_trim(all$detail)
all$desc   <- str_trim(all$desc)
saveRDS(all, "all.rds")

