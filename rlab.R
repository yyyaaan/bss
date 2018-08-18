
# simulate function -------------------------------------------------------

simData <- function(omega = matrix(c(5,6,7,7,9,8,1,2,3) , ncol = 3),
                    epsilon = matrix(c(1,2,3,3,2,5,6,1,3) * 1e-4, ncol = 3), 
                    N = 1e5,
                    model = list(ar = runif(2,-1,1))){
  
  p <- ncol(omega)
  if(p != ncol(epsilon)) stop("Epsilon and Omega must have same dimension")
  
  # simulate sources --------------------------------------------------------
  s <- replicate(3, arima.sim(model,N))
  s <- apply(s, 2, scale) # scaling to be zero-mean and unit-var: apply(s, 2, var); apply(s, 2, mean)

  # mixing ------------------------------------------------------------------
  x.simple <- s %*% t(omega)   # note for the matrix in transpose

  x.tv <- matrix(nrow = nrow(s), ncol = ncol(s))
  for (i in 1:N) 
    x.tv[i,] <- s[i,] %*% t(omega) %*% t(diag(3) + i * t(epsilon))
  
  # in case error, retry or re-model ----------------------------------------
  
  return(list(source = s, x.simple = x.simple, x.tv = x.tv, s.corr = cor(s)))
}

# sim <- simData(omega, epsilon, N, list(ar=runif(2,-1,1)))

# my JD -------------------------------------------------------------------

myCovJD <- function(covs){
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


genJD <- function(covs, method = "frjd"){
  # param covs as in JADE | in the form of dim = p,p,lag.max+1 | lag 0 at index 1
  # methods are frjd, fjd, rjd.fortran, FG, djd
  # NOTE for acf(.)$acf | dim = lag.max+1,p,p 
  if(dim(covs)[1] != dim(covs)[2]) stop('Autocovariance matrices should be in dim(p,p,k)')
  
  require(JADE)
  
  # whitening using S_0^{-1/2}, through eigen calculation
  eig <- eigen(covs[,,1])
  white <- eig$vectors %*% diag(eig$values^{-1/2}) %*% solve(eig$vectors)
  
  # whiten and format to rjd
  covs.white <- array(dim = c(dim(white),lag.max))
  for (i in 1:lag.max) covs.white[,,i] <- white %*% covs[,,i+1] %*% t(white)
  
  # joint diagnolization for V | note for the transpose of V
  # NOTE djd return type is different
  jd <- do.call(method, list(X = covs.white))
  if(method == "djd") {
    W <- jd %*% white
    D <- NA
  } else {
    W <- t(jd$V) %*% white
    D <- jd$D
  } 

  # D is the estimated diagnals
  list(W = W, D = D, white = white)
}



# sim by hand -------------------------------------------------------------

N <- 1e5
lag.max <- 4
sim <- list()

omega <- matrix(rnorm(9) , ncol = 3)
epsilon <- matrix(rnorm(9) * 1e-4, ncol = 3)
s1 <- arima.sim(list(ar=c(0.3,0.6)),N)
s2 <- arima.sim(list(ma=c(-0.3,0.3)),N)
s3 <- arima.sim(list(ar=c(-0.8,0.1)),N)
s <- cbind(s1,s2,s3)
s <- apply(s, 2, scale)

sim$x.simple <- s %*% t(omega)
sim$x.tv <- matrix(nrow = nrow(s), ncol = ncol(s))
for (i in 1:N) sim$x.tv[i,] <- s[i,] %*% t(omega) %*% t(diag(3) + i * t(epsilon))

      ### cleaning
remove(s1,s2,s3,s,i)









# test with original SOBI -------------------------------------------------

X <- sim$x.simple

  # autocovariance from acf
covs1 <- acf(X - colMeans(X), lag.max = lag.max, type = "covariance", plot = F)$acf
  # correct dimension to JADE
covs <- array(dim = dim(covs1)[3:1])
for (i in 0:lag.max) covs[,,i + 1] <- covs1[i + 1, ,]

MD(SOBI(X)$W, omega) # compare | SOBI usually better? why?
MD(myCovJD(covs)$W, omega)

      ### cleaning
remove(covs, covs1)

# TV-SOBI -----------------------------------------------------------------

# data("CPPdata"); X <- CPPdata
X <- sim$x.tv


lag.max <- 6
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

      ### clearning
remove(i,j,lag, x.design, y.design, est)

# step 2. JD for omega using Ra -------------------------------------------

jd <- myCovJD(Ra, method = "rjd")
W.est <- jd$W
omega.est <- solve(W.est)
MD(W.est, omega); MD(SOBI(X)$W, omega) # compared with normal SOBI


# step 3. solve for epsilon using Rb --------------------------------------


# Q in method 1 directly equals to Ra
Q <- Ra
# Q in method 2: Omega Lambda Omega
Q <- array(dim = c(p, p, lag.max + 1))
Q[,,1] <- omega.est %*% t(omega.est)
for (lag in 1:lag.max) Q[,,lag + 1] <- omega.est %*% jd$D[,,lag] %*%  t(omega.est)


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
for (i in 2:(lag.max + 1)) {
  y.design <- y.design + matrix(as.vector(Rb[,,i]), nrow = p^2)
  x.design <- x.design + H1[,,i] + H2[,,i]
}

est <- lm(y.design ~ x.design - 1)$coefficients
epsilon.est1 <- matrix(est, nrow = p, byrow = T)

      ### clearning
remove(Q, Qi, H1, H2, y.design, x.design, i, lag, est, pos, Ra, Rb, jd)


# step 3 alt: jd with Rc --------------------------------------------------

jd <- myCovJD(Rc)
  # thisW = (Epsilon * Omega)^{-1} = W * Epsion^{-1}
epsilon.est2 <- solve(jd$W) %*% W.est

      ### cleaning
remove(jd)

# results -----------------------------------------------------------------

  # tv-sobi
MD(W.est, omega)
MD(SOBI(X)$W, omega)
epsilon.est1 - epsilon; epsilon.est2 - epsilon

  # original sobi compare
MD(ans, omega)
MD(SOBI(sim$x.simple)$W, omega)

cat("Accuracy of Omega:\n using TV-SOBI", MD(W.est, omega), "\n using SOBI", MD(SOBI(X)$W, omega))


