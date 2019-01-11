

# TV-SOBI support ---------------------------------------------------------


nearestSPD <- function(X){
  if(length(dim(X)) != 2 ){
    stop("only works for matrix")
  } else if(dim(X)[1] != dim(X)[2]){
    stop("only works for square matrix")
  }
  
  RES <- Matrix::nearPD(X, doSym = TRUE)
  
  list(mat = RES$mat, normF = RES$normF)
}


approxJD <- function(covs, method = "frjd"){
  # follow the name in Miettinen, Nordhausen & Taskinen (2017)
  # param covs as in JADE | in the form of dim = p,p,lag.max+1 | lag 0 at index 1
  # NOTE for acf(.)$acf | dim = lag.max+1,p,p 
  
  require(JADE)
  
  if(dim(covs)[1] != dim(covs)[2]) 
    stop('Autocovariance matrices should be in dim(p,p,k+1)')
  
  lag.max <- dim(covs)[3] - 1
  

  # whitening using var^{-1/2}, through eigen calculation -------------------
  nearestDist <- 0
  EVD <- eigen(covs[ , , 1], symmetric = T)
  if (min(EVD$values) < 0) {
    # not semi-definite?
    spd <- nearestSPD(covs[ , , 1])
    EVD <- eigen(spd$mat, symmetric = T)
    warning(paste("covariance is NOT semi-definite, applying Nearest SPD;",
                  "Frobenius norm:", spd$normF))
    nearestDist <- spd$normF
  }
  white <- EVD$vectors %*% tcrossprod(diag(EVD$values^(-0.5)), EVD$vectors)
  

  # whiten autocovariance, format to rjd with symmetry fix ------------------
  covs.white <- array(dim = c(dim(white),lag.max))
  for (i in 1:lag.max) {
    covs.white[ , , i] <- white %*% tcrossprod(covs[ , , i+1], white)
    covs.white[ , , i] <- (covs.white[ , , i] + t(covs.white[ , , i]))/2
  }

  # joint diagnolization for V | note for the transpose of V
  JD <- do.call(method, list(X = covs.white))
  W <- crossprod(JD$V, white)
  W <- sweep(W, 1, sign(rowMeans(W)), "*")
  D <- JD$D

  # D is the estimated diagnals
  list(W = W, D = D, nearestDist = nearestDist)
}


# TV-SOBI Utilities -------------------------------------------------------

make_tvmix <- function(s, omega, epsilon){
  N <- nrow(s); p <- ncol(s)
  X <- matrix(nrow = N, ncol = p)
  for (i in 1:N){
    # X[i,] <- s[i,] %*% t(omega) %*% t(diag(2) + i * t(epsilon)) # actually equiv.
    X[i,] <- t( (diag(p) + i * epsilon) %*% omega %*% as.matrix(s[i,]) )
  } 
  return(X)
}

getW_t <- function(t, W, Epsilon){
  omega <- solve(W)
  p <- dim(omega)[1]
  if(is.null(Epsilon)) Epsilon <- diag(0, nrow = p) # compatible with SOBI
  if(is.na(as.vector(Epsilon)[1])) Epsilon <- diag(0, nrow = p) # compatible with SOBI
  omega_t <- (diag(p) + t * Epsilon) %*% omega
  return(solve(omega_t))
}

getMD_t <- function(omega_true, epsilon_true, t, W, Epsilon = NA){
  MD(getW_t(t, W, Epsilon), ((diag(dim(omega_true)[1]) + t * epsilon_true) %*% omega_true))
}

getMD_ave <- function(omega_true, epsilon_true, N, W, Epsilon = NA, ave_fun = function(x) mean(x)){
  mds <- numeric(N)
  for (i in 1:N) mds[i] <- getMD_t(omega_true, epsilon_true, i, W, Epsilon)
  return(ave_fun(mds))
}

restore_source <- function(X, W, Epsilon = NA){
  n <- nrow(X)
  z <- matrix(nrow = nrow(X), ncol = ncol(X))
  for (i in 1:n) z[i, ] <- X[i, ] %*% t(getW_t(i, W, Epsilon))
  return(z)
}

# TV-SOBI -----------------------------------------------------------------

tvsobi <- function(X, lag.max = 12, 
                   useQuadratic = TRUE,
                   epsilon.method = 1,
                   getSource = TRUE,
                   jd.method = "frjd"){

  if((!useQuadratic) & (epsilon.method == 3)) 
    stop("non-quadratic and the 3rd method for epsilon are incompatible")
  p <- ncol(X)
  n <- nrow(X)
  
  # step 1. estimate R0 R1 (R2) matrices ------------------------------------
  
  Ra <- Rb <- Rc <- array(dim = c(p, p, lag.max + 1))
  for (lag in 0:lag.max) {
    for (i in 1:p){
      for(j in 1:p){
        y.design <- X[(1 + lag):n, i] * X[1:(n - lag), j] # empirical autocovariance
        if( useQuadratic) h.design <- cbind(1:(n - lag), (1:(n - lag))^2)
        if(!useQuadratic) h.design <- 1:(n - lag)
        est <- lm(y.design ~ h.design)$coefficients # using lm
        Ra[i, j, lag + 1] <- est[1]
        Rb[i, j, lag + 1] <- est[2]
        if(useQuadratic) Rc[i, j, lag + 1] <- est[3]
      }
    }
  }
  # fix for symmetric
  for (i in 1:(lag.max + 1)) {
    Ra[ , , i] <- (Ra[ , , i] + t(Ra[ , , i]))/2
    Rb[ , , i] <- (Rb[ , , i] + t(Rb[ , , i]))/2
    if(useQuadratic) Rc[ , , i] <- (Rc[ , , i] + t(Rc[ , , i]))/2
  }
  remove(i, j, lag, h.design, y.design, est)
  
  
  
  # step 2. JD for omega using Ra -------------------------------------------
  
  jd <- approxJD(Ra, method = jd.method)
  nearestDist <- jd$nearestDist
  W.est <- jd$W
  omega.est <- solve(W.est)

  
  # step 3. solve for epsilon using Rb --------------------------------------

  if(epsilon.method %in% c(1,2)){
    if(epsilon.method == 2) { # = Omega Lambda Omega from JD in step 2
      Q <- array(dim = c(p, p, lag.max + 1))
      Q[,,1] <- omega.est %*% t(omega.est) 
      for (lag in 1:lag.max) Q[,,lag + 1] <- omega.est %*% jd$D[,,lag] %*%  t(omega.est)
    }
    
    if(epsilon.method == 1) Q <- Ra
    
    # design matricies
    H1 <- H2 <- array(0, dim = c(p^2, p^2, lag.max + 1))
    for(lag in 0:lag.max){
      for(i in 1:p^2){
        # the column to use  
        Qi <- Q[ , ceiling(i/p), lag+1]
        # H1 similar to diagnal
        pos <- (ifelse(i%%p == 0, p, i%%p ) - 1) * p  + (1:p)
        H1[i, pos, lag+1] <- Qi
        # H2 similar to vec
        pos <- (ceiling(i/p) - 1) * p + (1:p)
        H2[i, pos, lag+1] <- Qi
      }
    }
    
    y.design <-matrix(as.vector(Rb[ , , 1]), nrow = p^2)
    h.design <- H1[ , , 1] + H2[ , , 1]
    for (i in 2:(lag.max + 1)) {
      y.design <- y.design + matrix(as.vector(Rb[ , , i]), nrow = p^2)
      h.design <- h.design + H1[ , , i] + H2[ , , i]
    }
    
    est <- lm(y.design ~ h.design - 1)$coefficients
    epsilon.est <- matrix(est, nrow = p, byrow = T)
  } # end for method 1 and 2

  
  if(epsilon.method == 3) { # = jd(Rc)
    # this method directly return the results
    jd <- approxJD(Rc)
    epsilon.est <- NA
    epsilon.est <- try(solve(jd$W) %*% W.est, silent = T)
  }
  
  if(epsilon.method == 4) { # = jd(Rc)
    # this method directly return the results
    jd <- approxJD(Rc)
    epsilon.est <- NA
    epsilon.est <- try(W.est %*% solve(jd$W), silent = T)
  }
  
  
  if(!getSource)
    return(list(W = W.est, Epsilon = epsilon.est, nearestDist = nearestDist))
  
  S <- restore_source(X, W.est, epsilon.est)
  return(list(W = W.est, Epsilon = epsilon.est, S = S, nearestDist = nearestDist))
  
}