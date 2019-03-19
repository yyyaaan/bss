require(tidyverse)
require(JADE)
require(visNetwork)
require(RColorBrewer)

RColorBrewer::brewer.pal(3, "Set2")

# alogrithm alts ----------------------------------------------------------

processes <- c("Pre-processing and centering", "Decomposition of Autocovariances",
               "Yeredor TVSOBI", "LTVSOBI-1", "LTVSOBI-2", 
               "Restoration", "approxJD", "", "")
colors   <- RColorBrewer::brewer.pal(length(proc), "Set2")
outcomes  <- c("Observation", "Autocovariance Matrices", 
               "Beta_1", "Beta_2" , "Beta_3",
               "Restored Siganls",
               "Omega\n(Y-TVSOBI)", "Epsilon\n(Y-TVSOBI)",
               "Omega\n(LTVSOBI-1)", "Epsilon\n(LTVSOBI-1)",
               "Omega*Epsilon\n(LTVSOBI-2)","Omega\n(LTVSOBI-2)", "Epsilon\n(LTVSOBI-2)")
               
label_ids <- c(1, 9, 2, 9,
               7, 9, 9, 
               7, 9, 9, 9, 
               7, 9, 9, 9, 9, 
               9, 6, 9, 6, 9, 6)
edges <- data.frame(from = c(1, 2, 2, 2, 
                             7, 8, 8,
                             9, 9, 10, 10,
                             11, 12, 12, 13, 13,
                             7, 8, 9, 10, 12, 13),
                    to   = c(2, 3, 4 ,5, 
                             3, 4, 7,
                             3, 4, 4, 9,
                             5, 4, 11, 4, 11,
                             rep(6, 6)),
                    label = processes[label_ids], arrows = c(rep("to", 4), rep("from", 12), rep("to", 6)))
nodes <- data.frame(id = 1:length(outcomes), shape = "box",
                    label = outcomes, 
                    group = c(1, 1, 3, 3, 3, 1, 5, 5, 6, 6, 7, 7, 7),
                    level = c(1, 2, 3, 3, 3, 6, 4, 4.4, 4, 4.4, 4 ,4.4, 4.4))

visNetwork(nodes, edges) %>% visConfigure(enabled = T) %>% visHierarchicalLayout(direction = "UD", levelSeparation = 130) 



# performance -------------------------------------------------------------

plotIt <- function(sim, need_boxplot = F){

    # transform data ----------------------------------------------------------
  
  sim_tr <- sim %>% 
    mutate(sim_fun = str_remove(sim_fun, "function... sim_Generate_"), 
           bss_fun = str_remove(bss_fun, "function... ")) %>% 
    gather("key", "value", md_initial:sir_db) %>%
    filter(key %in% c("md_initial", "md_AVE", "sir_db"))
  
  sim_summary <- sim_tr %>%
    group_by(N, sim_fun, bss_fun, key) %>%
    summarise_at("value", mean, na.rm = TRUE)
  
  
  
    # mean and boxplots -------------------------------------------------------
  
  
  p <- sim_summary %>%
    filter(bss_fun %in% unique(.$bss_fun)) %>%
    ggplot(aes(x = N, y = value, color = bss_fun)) +
    geom_line() + geom_point(aes(shape = bss_fun)) +
    scale_x_log10(limits = c(min(sim$N), max(sim$N))) +
    facet_grid(key ~ sim_fun, scales = "free") +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = "Performance matrix (SignalType x Measurements)",
         subtitle = "Bootstrap-1000 result | color stands for different BSS method ")
  
  print(p)
  
  if(need_boxplot) {
    sim_tr %>%
      filter(bss_fun %in% unique(.$bss_fun)[3], sim_fun %in% unique(.$sim_fun)[1:2]) %>%
      ggplot(aes(x = as.factor(N), y = value, color = bss_fun)) +
      geom_boxplot() +
      facet_grid(key ~ sim_fun, scales = "free") +
      theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = "Performance matrix (SignalType x Measurements)",
           subtitle = "Bootstrap-500 result | selected measurements ")    
  }
}    


# ACTION ------------------------------------------------------------------


require(tidyverse); require(odbc)

sqlconn <- dbConnect(odbc::odbc(), "Study Database")
dbListTables(sqlconn)


sim <- sqlconn %>% dbGetQuery(paste0("select * from `", dbListTables(.)[3], "`"))

plotIt(sim)
