
# Nearest SPD -------------------------------------------------------------

nearestSPD <- function(X){
  if(length(dim(X)) != 2 ){
    stop("only works for matrix")
  } else if(dim(X)[1] != dim(X)[2]){
    stop("only works for square matrix")
  }
  
  RES <- Matrix::nearPD(X, ensureSymmetry = T)
  
  list(mat = RES$mat, normF = RES$normF)
}

# Approximate Joint Diagnolization ----------------------------------------
# follow the name in Miettinen, Nordhausen & Taskinen (2017)

approxJD <- function(covs, method = "frjd"){
  # param covs as in JADE | in the form of dim = p,p,lag.max+1 | lag 0 at index 1
  # NOTE for acf(.)$acf | dim = lag.max+1,p,p 
  
  require(JADE)
  
  if(dim(covs)[1] != dim(covs)[2]) 
    stop('Autocovariance matrices should be in dim(p,p,k+1)')
  lag.max <- dim(covs)[3] - 1
  
  # whitening using R_0^{-1/2}, through eigen calculation
  EVD <- eigen(covs[ , , 1], symmetric = T)
  if (min(EVD$values) < 0) {
    spd <- nearestSPD(covs[ , , 1])
    EVD <- eigen(spd$mat, symmetric = T)
    warning(paste("covariance is NOT semi-definite, applying Nearest SPD;",
                  "Frobenius norm:", spd$normF))
  }
  white <- EVD$vectors %*% tcrossprod(diag(EVD$values^(-0.5)), EVD$vectors)
  
  # whiten autocovariance, format to rjd with symmetry fix
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
  list(W = W, D = D)
}


# TV-SOBI -----------------------------------------------------------------

tvsobi <- function(X, lag.max = 12, 
                   useQuadratic = TRUE,
                   epsilon.method = 1,
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
  W.est <- jd$W
  omega.est <- solve(W.est)

  
  # step 3. solve for epsilon using Rb --------------------------------------

  if(epsilon.method == 3) { # = jd(Rc) and exit
    # this method directly return the results
    jd <- approxJD(Rc)
    epsilon.est <- NA
    epsilon.est <- try(solve(jd$W) %*% W.est, silent = T)
    return(list(W = W.est, Epsilon = epsilon.est, Ra = Ra, Rb = Rb, Rc = ifelse(is.na(Rc), NA, Rc)))
  }
  
  # for method 1 and 2, code continues after if{}
  if(epsilon.method == 2) { # = Omega Lambda Omega from JD in step 2
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
  # remove(Q, Qi, H1, H2, y.design, h.design, i, lag, est, pos, jd)
  
  return(list(W = W.est, Epsilon = epsilon.est, Ra = Ra, Rb = Rb, Rc = ifelse(is.na(Rc), NA, Rc)))
}



# printing detailed simulation --------------------------------------------



toAvoidRun <- function(){
  require(JADE)
  lag.max = 6

  # sim data ----------------------------------------------------------------
  N <- 1e5
  omega <- matrix(rnorm(9) , ncol = 3)
  epsilon <- matrix(rnorm(9) * 1e-7, ncol = 3)
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
  
  cat("\n\non simulated TIME-VARING mixture, MDIs:",
      "\nOriginal JADE::SOBI:\t",  MD(SOBI(X, lag.max)$W, omega),
      "\nSOBI with approxJD:\t",   MD(approxJD(covs)$W, omega),
      "\nNew TV-SOBI Quad-1:\t", MD(tvsobi(X, lag.max, useQuadratic = T, epsilon.method = 1)$W, omega),
      "\nNew TV-SOBI Quad-2:\t", MD(tvsobi(X, lag.max, useQuadratic = T, epsilon.method = 2)$W, omega),
      "\nNew TV-SOBI Quad-3:\t", MD(tvsobi(X, lag.max, useQuadratic = T, epsilon.method = 3)$W, omega),
      "\nNew TV-SOBI Line-1:\t", MD(tvsobi(X, lag.max, useQuadratic = F, epsilon.method = 1)$W, omega),
      "\nNew TV-SOBI Line-2:\t", MD(tvsobi(X, lag.max, useQuadratic = F, epsilon.method = 2)$W, omega),
      "\n")

  
  # original SOBI -----------------------------------------------------------
  
  # a normal mixture; to understand how approxJD works
  X <- z %*% t(omega)
  
  # autocovariance from acf, correct dimension to JADE
  covs1 <- acf(X - colMeans(X), lag.max = lag.max, type = "covariance", plot = F)$acf
  covs <- array(dim = dim(covs1)[3:1])
  for (i in 0:lag.max) covs[,,i + 1] <- covs1[i + 1, ,]
  
  cat("\n\non simulated ORDINARY mixture, MDIs:",
      "\nOriginal JADE::SOBI:\t",  MD(SOBI(X, lag.max)$W, omega),
      "\nSOBI with approxJD:\t",   MD(approxJD(covs)$W, omega),
      "\nNew TV-SOBI Quad-1:\t", MD(tvsobi(X, lag.max, useQuadratic = T, epsilon.method = 1)$W, omega),
      "\nNew TV-SOBI Quad-2:\t", MD(tvsobi(X, lag.max, useQuadratic = T, epsilon.method = 2)$W, omega),
      "\nNew TV-SOBI Quad-3:\t", MD(tvsobi(X, lag.max, useQuadratic = T, epsilon.method = 3)$W, omega),
      "\nNew TV-SOBI Line-1:\t", MD(tvsobi(X, lag.max, useQuadratic = F, epsilon.method = 1)$W, omega),
      "\nNew TV-SOBI Line-2:\t", MD(tvsobi(X, lag.max, useQuadratic = F, epsilon.method = 2)$W, omega),
      "\n")
}






# batch simulation recording ----------------------------------------------


yeredorSim <- function(N){
  
  require(JADE)

  # param
  omega <- matrix(c(3, -2, 1, 4), ncol = 2)
  epsilon <- matrix(c(-1, -2, 0.5, 1), ncol = 2) * 1e-4
  z <- rnorm(N)
  z1 <- 1 + 2*z^(-1) - 0.5*z^(-2) - z^(-3) + z^(-4)
  z <- rnorm(N)
  z2 <- 1 - z^(-1) + 3*z^(-2) + 2*z^(-3)
  z <- cbind(z1,z2)
#  z <- apply(cbind(z1,z2), 2, scale)
  X <- matrix(nrow = nrow(z), ncol = ncol(z))
  for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(2) + i * t(epsilon))

  start_time <- Sys.time()
  md_tvsobi <- tryCatch(MD(tvsobi(X)$W, omega), error = function(e) NA)
  time_tvsobi <- as.numeric(Sys.time() - start_time)
  
  start_time <- Sys.time() 
  md_sobi <- tryCatch(MD(SOBI(X)$W, omega), error = function(e) NA)
  time_sobi <- as.numeric(Sys.time() - start_time)
  
  
  c(md_tvsobi = md_tvsobi, md_sobi = md_sobi,
    time_tvsobi = time_tvsobi, time_sobi = time_sobi)
}


mySim <- function(N){
  require(JADE)
  lag.max = 6
  
  omega <- matrix(rnorm(9) , ncol = 3)
  epsilon <- matrix(rnorm(9) * 1e-5, ncol = 3)
  z1 <- arima.sim(list(ar=c(0.3,0.6)),N)
  z2 <- arima.sim(list(ma=c(-0.3,0.3)),N)
  z3 <- arima.sim(list(ar=c(-0.8,0.1)),N)
  z <- apply(cbind(z1,z2,z3), 2, scale)
  X <- matrix(nrow = nrow(z), ncol = ncol(z))
  for (i in 1:N) X[i,] <- z[i,] %*% t(omega) %*% t(diag(3) + i * t(epsilon))

  start_time <- Sys.time()
  md_tvsobi <- tryCatch(MD(tvsobi(X)$W, omega), error = function(e) NA)
  time_tvsobi <- as.numeric(Sys.time() - start_time)
  
  start_time <- Sys.time() 
  md_sobi <- tryCatch(MD(SOBI(X)$W, omega), error = function(e) NA)
  time_sobi <- as.numeric(Sys.time() - start_time)
  
  
  c(md_tvsobi = md_tvsobi, md_sobi = md_sobi,
    time_tvsobi = time_tvsobi, time_sobi = time_sobi)
}


simNplot <- function(){
  simResult1 <- data.frame(n = (1:100) * 100,
                           md_tvsobi = NA, md_sobi = NA,
                           time_tvsobi = NA, time_sobi = NA)
  simResult2 <- data.frame(n = (1:100) * 100,
                           md_tvsobi = NA, md_sobi = NA,
                           time_tvsobi = NA, time_sobi = NA)
  
  for(i in 1:nrow(simResult1)) {
    simResult1[i, 2:5] <- yeredorSim(simResult1$n[i])
    simResult2[i, 2:5] <- mySim(simResult2$n[i])
  }
  
  library(ggplot2)
  ggplot(simResult2) +
    geom_line(aes(n, md_tvsobi), color = "blue") + 
    geom_line(aes(n, md_sobi), color = "red")
}