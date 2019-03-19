require(tidyverse)
require(JADE)


res <- SOBI(fig_mixing$mix)
plot.ts(res$S)
plot.ts(fig_mixing$source)

# permutation -------------------------------------------------------------



# illustration of mixture -------------------------------------------------
N <- 1e4
omega <- matrix(rnorm(9) , ncol = 3)
epsilon <- matrix(rnorm(9) * 1e-4, ncol = 3)
z1 <- arima.sim(list(ar=c(0.3,0.6)),N)
z2 <- arima.sim(list(ma=c(-0.3,0.3)),N)
z3 <- arima.sim(list(ar=c(-0.8,0.1)),N)
z <- apply(cbind(z1,z2,z3), 2, scale)
X <- matrix(nrow = nrow(z), ncol = ncol(z))
for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(3) + i * t(epsilon))

  # z is normal mix, X is tv-mix
temp1 <- as.data.frame(X)
temp1$t <- 1:nrow(temp1)
temp1 <- melt(temp1, id.vars = "t")
temp1$type <- "time-varing mix"

temp2 <- as.data.frame(z)
colnames(temp2) <- c("V1", "V2", "V3")
temp2$t <- 1:nrow(temp2)
temp2 <- melt(temp2, id.vars = "t")
temp2$type <- "ordinary mix"

mix <- rbind(temp1, temp2)

ggplot(mix, aes(t, value, color = variable)) +
  geom_line() +
  facet_grid(variable~type) +
  theme_light()


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
