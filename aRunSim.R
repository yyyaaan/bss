remove(list = ls()) # clear environment
source("rfun.R"); source("rsim.R"); source("rlab_x.R");
library(tidyverse)


# simulation sampling interval --------------------------------------------


do_it_once <- function(x, z, lll = 6, id = "hi"){
  
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


path <- paste0(getwd(), "/sim/")

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
  if (file.exists(fname)) {
    saveRDS(rbind(readRDS(fname), df), file = fname)
    
  } else {
    saveRDS(df, file = fname)
  }
}


# call --------------------------------------------------------------------

# p <- 3
# Omega   <- matrix(runif(p^2, -10, 10), ncol = p)
# Epsilon <-  1e-3 * matrix(runif(p^2, -10, 10), ncol = p)

Omega <- matrix(c(2, -9, -4, -6, 5, 6, 0.5, 3, 8), ncol =3)
Epsilon <- 1e-5 * matrix(c(-3, -4, 9, 6, 2.5, 2.1, -6, 6, 7), ncol = 3)
zall <- sim_good_sources(N = 1e5, 3)
xall <- tvmix(zall, Omega, Epsilon)
save(Omega, Epsilon, xall, zall, file = paste0(getwd(),"/sim/E5N5_A.rdata"))

# loop for freqs
freq_list <- 2^(0:10)
for(freq in freq_list){
  for(l in c(3,6,12,1)){
    ids  <- seq(from = 1, to = nrow(xall), by = freq)
    x <- xall[ids,]
    z <- zall[ids,]
    do_it_once(x, z, lll = l, id = paste0("fixed_freq_E5N5_A_lag", l))
  }
}




# boosting ----------------------------------------------------------------


for(i in 1:100){
  Omega <- matrix(c(2, -9, -4, -6, 5, 6, 0.5, 3, 8), ncol =3)
  Epsilon <- 1e-5 * matrix(c(-3, -4, 9, 6, 2.5, 2.1, -6, 6, 7), ncol = 3)
  zall <- sim_good_sources(N = 1e5, 3)
  xall <- tvmix(zall, Omega, Epsilon)
  
  # loop for freqs
  freq_list <- 2^(0:10)
  for(freq in freq_list){
    for(l in c(3,6,12,1)){
      ids  <- seq(from = 1, to = nrow(xall), by = freq)
      x <- xall[ids,]
      z <- zall[ids,]
      do_it_once(x, z, lll = l, id = paste0("fixed_freq_E5N5_Boot_lag", l))
    }
  }
}




# results of boosting -----------------------------------------------------

library(tidyverse)
res <- readRDS("/home/yan/bss/aBootRes.rds")

res_sum <- res %>% 
  mutate(series = str_sub(id, 12, 15)) %>%
  group_by(criteria, series, method, N, p, lag) %>%
  summarise_at("value", mean) 

for (i in 1:nrow(res_sum)){
if(res_sum$series[i] == "E4N4") res_sum$seriesT[i] = "Set IV"
if(res_sum$series[i] == "E4N5") res_sum$seriesT[i] = "Set II"
if(res_sum$series[i] == "E5N4") res_sum$seriesT[i] = "Set III"
if(res_sum$series[i] == "E5N5") res_sum$seriesT[i] = "Set I"

if(res_sum$lag[i] == 1) res_sum$lagT[i] = "Lag = 1"
if(res_sum$lag[i] == 3) res_sum$lagT[i] = "Lag = 3"
if(res_sum$lag[i] == 6) res_sum$lagT[i] = "Lag = 6"
if(res_sum$lag[i] == 12) res_sum$lagT[i] = "Lag = 12"
}

saveRDS(res_sum, file = "/home/yan/bss/thesis/bss_res.rds")

m <- unique(res_sum$criteria)[4]

readRDS("/home/yan/bss/thesis/bss_res.rds") %>%
  filter(criteria == m, N > 49, method != "frjd", lag != 1, seriesT %in% c("Set I", "Set II")) %>%
  ggplot(aes(N, value, color = method)) +
  geom_point() + geom_line() + 
  scale_x_log10() +
  facet_grid(seriesT~lagT) +
  theme_light()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


