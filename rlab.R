
N <- 1e3
lag.max <- 12

# simulate function -------------------------------------------------------

simData <- function(omega = matrix(c(5,6,7,7,9,8,1,2,3) , ncol = 3),
                    epsilon = matrix(c(1,2,3,3,2,5,6,1,3) * 1e-4, ncol = 3),
                    N = 1e3, plot = FALSE, as.ts = FALSE){
  
  p <- ncol(omega)
  if(p != ncol(epsilon)) stop("Epsilon and Omega must have same dimension")
  
  # simulate sources --------------------------------------------------------
  s <- matrix(nrow = N, ncol = p)
  for (i in 1:p){
    s.i <- arima.sim(list(ar = runif(1)), N)
    s[,i] <- scale(s.i)
  }  
  # scaling to be zero-mean and unit-var: apply(s, 2, var); apply(s, 2, mean)
  if(plot) plot.ts(s)
  
  # mixing ------------------------------------------------------------------
  
  x.simple <- s %*% t(omega)   # note for the matrix in transpose
  if(plot) plot.ts(x.simple)
  
  x.tv <- matrix(nrow = nrow(s), ncol = ncol(s))
  for (i in 1:N) 
    x.tv[i,] <- s[i,] %*% t(omega) %*% t(diag(3) + i * t(epsilon))
  if(plot) plot.ts(x.tv)
  
  return(list(source = s, x.simple = x.simple, x.tv = x.tv))
}



# simulate data -----------------------------------------------------------

epsilon <- matrix(c(1,2,3,3,2,5,6,1,3) * 1e-4, ncol = 3)
omega <- matrix(c(5,6,7,7,9,8,1,2,3) , ncol = 3)
sim <- simData(omega, epsilon, plot = T)

# SOBI in detail ----------------------------------------------------------

myCovJD <- function(covs){
  # param covs as in acf(.)$acf | dim = p,p,lag.max+1 including lag=0
  require(JADE)
  
  # whitening using S_0^{-1/2}, through eigen calculation
  eig <- eigen(covs[1,,])
  white <- solve(eig$vectors %*% sqrt(diag(eig$values)))
  
  # whiten and format to rjd
  covs.white <- array(dim = c(dim(white),lag.max))
  for (i in 1:lag.max) covs.white[,,i] <- white %*% covs[i+1,,] %*% t(white)
  
  # joint diagnolization for V | note for the transpose of V
  V <- frjd(covs.white, maxiter = 1e9)$V
  t(V) %*% white
}

X <- sim$x.simple
ans1 <- myCovJD(acf(X, lag.max = lag.max, type = "covariance", plot = F)$acf)
ans0 <- SOBI(X)$W

# compare results | SOBI better? why?
MD(ans1, omega); MD(ans0, omega)
print(ans1); print(ans0)


# TV-SOBI -----------------------------------------------------------------


X <- sim$x.tv
p <- ncol(X)
n <- nrow(X)
R0 <- R1 <- R2 <- array(dim = c(lag.max + 1, p, p))

for (lag in 0:lag.max) {
  for (i in 1:p){
    for(j in 1:p){
      # lag = 1; i = 2; j = 3
      y.design <- X[(1+lag):n,i] * X[1:(n-lag),j] # empirical autocovariance
      x.design <- cbind(rep(1, n-lag), 1:(n-lag), (1:(n-lag))^2)
      est <- lm(y.design ~ x.design - 1)$coefficients # using lm avoids the singular problem
      R0[lag+1, i, j] <- est[1]
      R1[lag+1, i, j] <- est[2]
      R2[lag+1, i, j] <- est[3]
    }
  }
}

MD(myCovJD(R0), omega)
R2



