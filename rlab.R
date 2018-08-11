T <- 500
N <- 2*T + 1

# generate data -----------------------------------------------------------

df <- data.frame(z1 = rnorm(N), z2 = rt(N, 2), z3 = runif(N))
z.df <- scale(df, scale = F) 
  ### must be centered
z.ts <- ts(z.df, start =  - T)
z <- t(as.matrix(z.df)) 
  ### R data frame - it is transposed!
plot(z.ts)



# acf, stationary
cov.tau <- acf(z.ts, lag.max = 12, type = "covariance")$acf
cov.tau[1,,] 
  ### cov, at lag 0, there is 1 step difference


z[,2:1001] %*% t(z[,1:1000]) / N


a <- acf(z.ts, lag.max = 12)
a$lag
