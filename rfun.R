# rstudioapi::viewer("https://boring.fi/bss")
# x %*% t(y) (tcrossprod).
# follow [p,p,NN]

library(magrittr)

# utilities ---------------------------------------------------------------

tvmix <- function(z, Omega, Epsilon, x_only = T){
  N <- nrow(z);  p <- ncol(z)
  # mean correction
  zz <- scale(z, center = T, scale = F)
  # mixing with time index t
  x <- matrix(ncol = p, nrow = N)
  Omega_t <- array(dim = c(p,p,N))
  for (t in 1:N) {
    Omega_t[ , , t] <- (diag(p) + t * Epsilon) %*% Omega
    x[t, ] <- zz[t, ] %*% t(Omega_t[, , t])
  }
  # ready
  if(x_only) return(x)
  return(list(x = x, z = z, Omega = Omega, Epsilon = Epsilon, 
              Omega_t = Omega_t, means = colMeans(z)))
}

# pseudo data -------------------------------------------------------------

N <- 1e3
z <-cbind(var1 = arima.sim(list(ar=runif(1,-1,1)),N),
          var2 = arima.sim(list(ar=runif(1,-1,1)),N),
          var3 = arima.sim(list(ar=runif(1,-1,1)),N))
Omega <- matrix(runif(9, 1, 10), ncol = 3)
Epsilon <- 1e-4 * matrix(runif(9, 1, 10), ncol = 3)
lags <- 6
x <- tvmix(z, Omega, Epsilon)

lapply(step1(x, lags), mean)


# algorithm ---------------------------------------------------------------


# linear estimation for covariance separation
# i <- 2; j <- 2; l <- 5; lm_items = 3
step1 <- function(x, lags){
  N <- nrow(x); p <- ncol(x)

    # lags accept two types | we need 0 for whitening
  if(length(lags) == 1) lags <- 0:lags 
  if(!(0 %in% lags)) lags <- c(0, lags)
  
    # looping for each lags and element-wise matrix
  empty_mat <- array(dim = c(p,p, length(lags)))
  beta_1 <- beta_2 <- beta_3 <- r_squared <- rep(list(empty_mat), 5); # array(dim = c(p,p, length(lags)))
  
  for(l in lags) { for (i in 1:p) { for (j in 1:p) {

          # build y
        y <- numeric(N-l)
        for (t in 1:(N-l)){
          cov_x_t_l <- matrix(x[t, ], nrow = 3) %*% matrix(x[t + l, ], ncol = 3)
          y[t] <- cov_x_t_l[i, j]
        }
          # build H, except for 1 columns
        H2 <-  1:(N-l)
        H3a <- (1:(N-l))^2
        H3b <- (1:(N-l))^2 + l
        H3c <- (1:(N-l))^2 + (1:(N-l))
        H3d <- (1:(N-l))^2 + (1:(N-l)) * l
        lm_data = data.frame(y, H2, H3a, H3b, H3c, H3d)
        
          # regression for betas in alternative ways
        lm_res <- list()
        lm_res[[1]] <- lm(y ~ H2, lm_data)
        lm_res[[2]] <- lm(y ~ H2 + H3a, lm_data)
        lm_res[[3]] <- lm(y ~ H2 + H3b, lm_data)
        lm_res[[4]] <- lm(y ~ H2 + H3c, lm_data)
        lm_res[[5]] <- lm(y ~ H2 + H3d, lm_data)

          # get result
        for(alt in 1:6){
          r_squared[[alt]][i,j,l+1] <- summary(lm_res[[alt]])$r.squared
          beta_1[[alt]][i,j,l+1] <- lm_res[[alt]]$coefficients[1]
          beta_2[[alt]][i,j,l+1] <- lm_res[[alt]]$coefficients[2]
          beta_3[[alt]][i,j,l+1] <- lm_res[[alt]]$coefficients[3]
        }
    }}} # end loop for i,j,l
  
  r_squared
}




step2 <- function(covs, method = "frjd"){
  # ApproxJD, Miettinen, Nordhausen & Taskinen (2017)
  # param covs as in JADE | in the form of dim = p,p,lag.max+1 | lag 0 at index 1
  # NOTE for acf(.)$acf | dim = lag.max+1,p,p 
  
  require(JADE)
  
  if(dim(covs)[1] != dim(covs)[2]) 
    stop('Autocovariance matrices should be in dim(p,p,k+1)')

    # here, for sure 0 should be removed
  lag_len <- dim(covs)[3] - 1 
  
    # whitening using var^{-1/2}, through eigen calculation
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
  
    # whiten autocovariance, format to rjd with symmetry fix
  covs.white <- array(dim = c(dim(white),lag_len))
  for (i in 1:lag_len) {
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


r_squared %>% mean()
