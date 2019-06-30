remove(list = ls()) # clear environment
source("rfun.R"); source("rlab_x.R"); source("rsim.R")
library(tidyverse)
options(stringsAsFactors = FALSE)

# results of boosting -----------------------------------------------------
# this supports aRunSim.R

res <- readRDS("/home/yan/bss/aBootRes.rds")

res %>% 
  mutate(series = str_sub(id, 12, 15)) %>%
  filter(criteria == "SIR_diag_sq", N > 90, method == "LTV-SOBI", lag != 1, series == "E5N5") %>%
  ggplot(aes(as.factor(N), value)) + geom_boxplot(color = "darkgrey")  +
  scale_y_continuous(name="tvSIR") +
  scale_x_discrete(name="Sampling Size") +
  theme_light()


res_sum <- res %>% 
  mutate(series = str_sub(id, 12, 15)) %>%
  group_by(criteria, series, method, N, p, lag) %>%
  summarise_at("value", mean) 

res_N <- res %>% 
  mutate(series = str_sub(id, 12, 15)) %>%
  group_by(criteria, series, desc, N, p, lag) %>%
  count()

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

# old benchmarking plotting -----------------------------------------------

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




sum(duplicated(rbind(aBootRes, all)))
