N <- 1e5
lag.max <- 4

# simulate function -------------------------------------------------------

simData <- function(omega = matrix(c(5,6,7,7,9,8,1,2,3) , ncol = 3),
                    epsilon = diag(3), N = 1e3, plot = FALSE, as.ts = FALSE){
  
  p <- ncol(omega)
  if(p != ncol(epsilon)) stop("Epsilon and Omega must have same dimension")
  
  # simulate sources --------------------------------------------------------
  sim <- function() arima.sim(list(ar = runif(2,-1,1)), N)
  
  s <- matrix(nrow = N, ncol = p)
  for (i in 1:p){
    s.i <- tryCatch(sim(), error = function(e) sim())
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
  
  # because of stationality, it may return error. In such case, retry
  return(list(source = s, x.simple = x.simple, x.tv = x.tv))
}




# SOBI in detail ----------------------------------------------------------

myCovJD <- function(covs, style = ""){
  # param covs as in JADE | in the form of dim = p,p,lag.max+1 | lag 0 at index 1
  # NOTE for acf(.)$acf | dim = lag.max+1,p,p 
  if(dim(covs)[1] != dim(covs)[2]) stop('Autocovariance matrices should be in dim(p,p,k)')

  require(JADE)
  
  # whitening using S_0^{-1/2}, through eigen calculation
  eig <- eigen(covs[,,1])
  white <- solve(eig$vectors %*% sqrt(diag(eig$values)))
  
  # whiten and format to rjd
  covs.white <- array(dim = c(dim(white),lag.max))
  for (i in 1:lag.max) covs.white[,,i] <- white %*% covs[,,i+1] %*% t(white)
  
  # joint diagnolization for V | note for the transpose of V
  jd <- frjd(covs.white, maxiter = 1e5)
  
  # D is the estimated diagnals
  list(W = t(jd$V) %*% white, D = jd$D, white = white)
}



# test with original SOBI -------------------------------------------------


epsilon <- matrix(c(1,2,3,3,2,5,6,1,3) * 1e-4, ncol = 3)
omega <- matrix(c(5,6,7,7,9,8,1,2,3) , ncol = 3)
sim <- simData(omega, epsilon, plot = F)
X <- sim$x.simple
  
  # autocovariance from acf
covs1 <- acf(X, lag.max = lag.max, type = "covariance", plot = F)$acf
  # correct dimension to JADE
covs <- array(dim = dim(covs1)[3:1])
for (i in 0:lag.max) covs[,,i + 1] <- covs1[i + 1, ,]

ans <- myCovJD(covs)$W
round(ans %*% omega, 2) # roughly permutation?
MD(ans, omega); MD(SOBI(X)$W, omega) # compare | SOBI usually better? why?

# TV-SOBI -----------------------------------------------------------------


epsilon <- matrix(c(1,2,3,3,2,5,6,1,3) * 1e-4, ncol = 3)
omega <- matrix(c(5,6,7,7,9,8,1,2,3) , ncol = 3)
sim <- simData(omega, epsilon, plot = F)
X <- sim$x.tv

p <- ncol(X)
n <- nrow(X)


# step 1. estimate R0 R1 R2 matrices --------------------------------------

Ra <- Rb <- Rc <- array(dim = c(p, p, lag.max+1))

for (lag in 0:lag.max) {
  for (i in 1:p){
    for(j in 1:p){
      # lag = 1; i = 2; j = 3
      y.design <- X[(1+lag):n,i] * X[1:(n-lag),j] # empirical autocovariance
      x.design <- cbind(rep(1, n-lag), 1:(n-lag), (1:(n-lag))^2)
      est <- lm(y.design ~ x.design - 1)$coefficients # using lm avoids the singular problem
      Ra[i, j, lag+1] <- est[1]
      Rb[i, j, lag+1] <- est[2]
      Rc[i, j, lag+1] <- est[3]
    }
  }
}

  # fix for symmetric
for (i in 1:(lag.max+1)) {
  Ra[,,i] <- (Ra[,,i] + t(Ra[,,i]))/2
  Rb[,,i] <- (Rb[,,i] + t(Rb[,,i]))/2
  Rc[,,i] <- (Rc[,,i] + t(Rc[,,i]))/2
}


# step 2. JD for omega using Ra -------------------------------------------

jd <- myCovJD(Ra)
W.est <- jd$W
omega.est <- solve(W.est)
MD(W.est, omega); MD(SOBI(X)$W, omega) # compared with normal SOBI


# step 3. solve for epsilon using Rb --------------------------------------

  # Q in method 1 directly equals to Ra
Q1 <- Ra
  # Q in method 2: Omega Lambda Omega
Q2 <- array(dim = c(p, p, lag.max + 1))
Q2[,,1] <- omega.est %*% t(omega.est)
for (lag in 1:lag.max) {
  Q2[,,lag + 1] <- omega.est %*% jd$D[,,lag] %*%  t(omega.est)
} 
  # select Q | the difference is actually negligible, depending on JD result
Q <- Q1

  # design matricies
H1 <- H2 <- array(0, dim = c(p^2, p^2, lag.max + 1))
for(lag in 0:lag.max){
  for(i in 1:p^2){
      # the column to use  
    Qi <- Q[,ceiling(i/p),lag+1]
      # H1 similar to diagnal
    pos <- (ifelse(i%%p == 0, p, i%%p ) - 1) * p  + (1:p)
    H1[i, pos, lag+1] <- Qi
      # H2 similar to vec
    pos <- (ceiling(i/p) - 1) * p + (1:p)
    H2[i, pos, lag+1] <- Qi
  }
}

y.design <-matrix(as.vector(Rb[,,1]), nrow = p^2)
x.design <- H1[,,1] + H2[,,1]
for (i in 2:p) {
  y.design <- y.design + matrix(as.vector(Rb[,,i]), nrow = p^2)
  x.design <- x.design + H1[,,i] + H2[,,i]
}

est <- lm(y.design ~ x.design - 1)$coefficients
epsilon.est <- matrix(est, nrow = p)
MD(epsilon.est, epsilon)
