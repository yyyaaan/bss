# plotting only -----------------------------------------------------------
library(reshape2); library(ggplot2)

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

SimRes <- readRDS("zzz_sim_batch.rds")
SimRes$n <- as.numeric(SimRes$n)

SimRes %>% 
#  filter(flag == 0) %>%
  group_by(method, n, simtype) %>%
  summarise_at("md", mean) %>%
  ggplot(aes(x = n, y = md, group = method,color = method)) +
  geom_point() + geom_line() +
  scale_x_log10(limits = c(min(SimRes$n), max(SimRes$n))) +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Sample Size n", y = "Mean MD-value") +
  facet_grid(simtype ~ .) +
  theme_light()