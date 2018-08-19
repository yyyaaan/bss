# use git, especially v1.0 for older versions
# removed from 1.0 simData FUN and other simulation
# add nonJD FUN
# add tvsobi FUN
# simulation example add to function to avoid run

# non-orthogonal JD -------------------------------------------------------

nonJD <- function(covs, method = "frjd", fix = FALSE){
  # param covs as in JADE | in the form of dim = p,p,lag.max+1 | lag 0 at index 1
  # methods are frjd, fjd, djd
  # NOTE for acf(.)$acf | dim = lag.max+1,p,p 
  
  require(JADE)
  
  if(dim(covs)[1] != dim(covs)[2]) stop('Autocovariance matrices should be in dim(p,p,k)')
  lag.max <- dim(covs)[3] - 1
  
  # whitening using R_0^{-1/2}, through eigen calculation
  # error handling is needed
  eig <- eigen(covs[,,1])
  if (min(eig$values) < 0) {
    if(!fix) stop("first matrix used for whitening is NOT postive semi-definite\nAdd fix = TRUE can dismiss the error but at cost of accuracy")
    eig$values[eig$values<0] <- 0
  }
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

# tv-sobi -----------------------------------------------------------------

tvsobi <- function(X, lag.max = 12, jd.method = "frjd", epsilon.method = 1, fix = F){
  p <- ncol(X)
  n <- nrow(X)
  
  # step 1. estimate R0 R1 R2 matrices --------------------------------------
  
  Ra <- Rb <- Rc <- array(dim = c(p, p, lag.max+1))
  for (lag in 0:lag.max) {
    for (i in 1:p){
      for(j in 1:p){
        # lag = 1; i = 2; j = 3
        y.design <- X[(1+lag):n,i] * X[1:(n-lag),j] # empirical autocovariance
        h.design <- cbind(rep(1, n-lag), 1:(n-lag), (1:(n-lag))^2)
        est <- lm(y.design ~ h.design - 1)$coefficients # using lm avoids the singular problem
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
  # remove useless
  remove(i,j,lag, h.design, y.design, est)
  
  # step 2. JD for omega using Ra -------------------------------------------
  jd <- nonJD(Ra, method = jd.method, fix = fix)
  W.est <- jd$W
  omega.est <- solve(W.est)

  # step 3. solve for epsilon using Rb --------------------------------------
    
  if(epsilon.method == 3) { # = jd(Rc) and exit
    # this method directly return the results
    # potential error -> warnings and return NA for Epsilon
    jd <- tryCatch(nonJD(Rc), 
                   error = function(e) warning("Rc[,,1] is NOT postive semi-definite\nAdd fix = TRUE can dismiss the error but at cost of accuracy"))
    epsilon.est <- NA
    epsilon.est <- try(solve(jd$W) %*% W.est, silent = T)
    
    return(list(W = W.est, Epsilon = epsilon.est, Ra = Ra, Rb = Rb, Rc = Rc))
  }
  
  # for method 1 and 2, code continues after if{}
  if(epsilon.method == 2) { # = Omega Lambda Omega from JD in step 2
    if(jd.method == "djd") {
      warning("djd is not compatible with Epsilon estimation method 2")
      return(list(W = W.est, Epsilon = NA, Ra = Ra, Rb = Rb, Rc = Rc))
    }
    
    Q <- array(dim = c(p, p, lag.max + 1))
    Q[,,1] <- omega.est %*% t(omega.est) 
    for (lag in 1:lag.max) Q[,,lag + 1] <- omega.est %*% jd$D[,,lag] %*%  t(omega.est)
  }

  if(epsilon.method == 1) { # = Ra
    Q <- Ra
  }
  
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
  h.design <- H1[,,1] + H2[,,1]
  for (i in 2:(lag.max + 1)) {
    y.design <- y.design + matrix(as.vector(Rb[,,i]), nrow = p^2)
    h.design <- h.design + H1[,,i] + H2[,,i]
  }
  
  est <- lm(y.design ~ h.design - 1)$coefficients
  epsilon.est <- matrix(est, nrow = p, byrow = T)
  # remove(Q, Qi, H1, H2, y.design, h.design, i, lag, est, pos, jd)
  
  return(list(W = W.est, Epsilon = epsilon.est, Ra = Ra, Rb = Rb, Rc = Rc))
}



# simulation test ---------------------------------------------------------

toAvoidRun <- function(){
  
  lag.max = 6

  # sim data ----------------------------------------------------------------
  N <- 1e5
  omega <- matrix(rnorm(9) , ncol = 3)
  epsilon <- matrix(rnorm(9) * 1e-4, ncol = 3)
  z1 <- arima.sim(list(ar=c(0.3,0.6)),N)
  z2 <- arima.sim(list(ma=c(-0.3,0.3)),N)
  z3 <- arima.sim(list(ar=c(-0.8,0.1)),N)
  z <- apply(cbind(z1,z2,z3), 2, scale)
  X <- matrix(nrow = nrow(z), ncol = ncol(z))
  for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(3) + i * t(epsilon))


  # tv-SOBI -----------------------------------------------------------------
  
  covs1 <- acf(X, lag.max = lag.max, type = "covariance", plot = F)$acf
  covs <- array(dim = dim(covs1)[3:1])
  for (i in 0:lag.max) covs[,,i + 1] <- covs1[i + 1, ,]
  
  cat("\n\non simulated time varying mixture, the minimum distance index:",
      "\nOriginal SOBI:",  MD(SOBI(X, lag.max)$W, omega),
      "\nSOBI with nonJD:",   MD(nonJD(covs)$W, omega),
      "\nNew TV-SOBI:", MD(tvsobi(X, lag.max)$W, omega))

  
  # original SOBI -----------------------------------------------------------
  
  # a normal mixture; to understand how nonJD works
  X <- z %*% t(omega)
  
  # autocovariance from acf, correct dimension to JADE
  covs1 <- acf(X - colMeans(X), lag.max = lag.max, type = "covariance", plot = F)$acf
  covs <- array(dim = dim(covs1)[3:1])
  for (i in 0:lag.max) covs[,,i + 1] <- covs1[i + 1, ,]
  
  cat("\n\non simulated ordinary mixture, the minimum distance index:",
      "\nOriginal SOBI:",  MD(SOBI(X, lag.max)$W, omega),
      "\nSOBI with nonJD:",   MD(nonJD(covs)$W, omega),
      "\nNew TV-SOBI:", MD(tvsobi(X, lag.max)$W, omega))
}









