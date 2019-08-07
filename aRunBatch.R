remove(list = ls()) # clear environment
source("rfun.R"); source("rlab_x.R"); source("rsim.R")
library(tidyverse)
options(stringsAsFactors = FALSE)


# results of boosting -----------------------------------------------------
# this supports aRunSim.R
library(tidyverse)
aBootRes <- readRDS("/home/yanpan/bss/aBootRes.rds")
aBootClean <- aBootRes %>% 
  mutate(series  = case_when(substr(id, 7, 10) == "_fix" ~ substr(id, 19, 22),
                             substr(id, 5, 8)  == "_fix" ~ substr(id, 17, 20),
                             substr(id, 1, 4)  == "fixe" ~ substr(id, 12, 15))) %>%
  mutate(seriesT = case_when(series == "E4N4" ~ "Set IV",
                             series == "E4N5" ~ "Set II",
                             series == "E5N4" ~ "Set III",
                             series == "E5N5" ~ "Set I"),
         lagT    = case_when(lag == 1  ~ "Lag = 1",
                             lag == 3  ~ "Lag = 3",
                             lag == 6  ~ "Lag = 6",
                             lag == 12 ~ "Lag = 12"),
         detail  = str_replace_all(detail, "YeredorTVOBI", "YeredorTVSOBI"))
aBootClean <- aBootClean %>%
  mutate(desc   = str_replace_all(desc  , "YeredorTVOBI", "YeredorTVSOBI"),
         method = str_replace_all(method, "YeredorTVOBI", "YeredorTVSOBI"),
         lagT   = factor(lagT, levels = c("Lag = 1", "Lag = 3", "Lag = 6", "Lag = 12"))) %>%
  group_by(criteria, seriesT, method, N, p, lagT)

saveRDS(aBootClean, "aBootClean.rds")
res <- aBootClean
save(res, fig_mixing, file = "/home/yanpan/bss/thesis/thesis.rdata", compress = TRUE)
# res_sum2 <- aBootClean %>% summarise_at("value", mean)
# res_N    <- aBootClean %>% count()



# short handed for detailed plot ------------------------------------------
res_sum <- res %>% 
  mutate(series = str_sub(id, 12, 15)) %>%
  group_by(criteria, series, detail, N, p, lag) %>%
  summarise_at("value", mean) 
  # do the loop above
res_sum$method <- res_sum$detail
res_sum$flag <- substr(res_sum$method, 1, 9)

m <- unique(res_sum2$criteria)[4]; 
cutoff <- 90; setVector <- c("Set I", "Set II"); yaxisname <- "tvSIR"; limitsVector <- c(0,9);
res_sum %>%
  filter(flag == "LTV-SOBI ", criteria == m, N > cutoff, method != "frjd", lag != 1, seriesT %in% setVector) %>%
  ggplot(aes(N, value, color=method, shape = method, linetype = method)) +
  #geom_point(size = 1) + 
  geom_line() + 
  scale_color_brewer(palette = "Set2") +
  scale_x_log10() +
  facet_grid(seriesT~lagT) +
  scale_y_continuous(name=yaxisname, limits=limitsVector) +
  theme(legend.position="bottom",
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text.x = element_text(angle = 90))


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

