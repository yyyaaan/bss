require(tidyverse)

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

require(tidyverse)

sim <- readRDS("res_sim_boot.rds") %>% 
  filter(method != "TVSOBI, Quadratic FALSE , Epsilon Mehtod 3")



sim_summary <- sim %>%
  group_by(N, p, method) %>%
  summarise_at(c(4:which(colnames(sim) == "md_AVE")), function(x) mean(x, na.rm = TRUE)) %>% 
  gather(4:which(colnames(sim) == "md_AVE"), key = "measured_time", value = "MD")

  # plot for mean values ----------------------------------------------------

sim_summary %>%
  filter(p == 2) %>%
  ggplot(aes(x = N, y = MD, color = measured_time)) +
  geom_line() +
  scale_x_log10(limits = c(min(sim$N), max(sim$N))) +
  facet_wrap(~method) +
  ggtitle("MD measure at different time (bootstrap result)", subtitle = "Yeredor's simulation (2-dimensional)")
  
sim_summary %>%
  filter(p == 3) %>%
  ggplot(aes(x = N, y = MD, color = measured_time)) +
  geom_line() +
  scale_x_log10(limits = c(min(sim$N), max(sim$N))) +
  facet_wrap(~method) +
  ggtitle("MD measure at different time (bootstrap result)", subtitle = "ARIMA simulation (3-dimensional)")


  # boxplots ----------------------------------------------------------------

sim$N <- as.factor(sim$N)
sim$p <- as.factor(sim$p)

sim %>%
  ggplot(aes(x = N, y = md_AVE)) +
  geom_boxplot(aes(color = p)) +
  facet_grid(p~method) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Average MD over all time")

sim %>%
  ggplot(aes(x = N, y = md_initial)) +
  geom_boxplot(aes(color = p)) +
  facet_grid(p~method) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("MD at t=0 (i.e. Omega only)")


sim$epsilon[999] %>% str_split(",", simplify = T) %>% as.numeric() %>% matrix(ncol = sqrt(length(.)))
