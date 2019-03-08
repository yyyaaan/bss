# rstudioapi::viewer("https://boring.fi/bss")
# x %*% t(y) (tcrossprod).
# follow [p,p,NN]
## as.vector = vec()

library(magrittr)

# utilities ---------------------------------------------------------------

tvmix <- function(z, Omega, Epsilon, x_only = TRUE){
  N <- nrow(z);  p <- ncol(z)
  # mean correction
  zz <- scale(z, center = TRUE, scale = FALSE)
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

tvunmix <- function(x, Omega_hat, Epsilon_hat){
  N <- nrow(z);  p <- ncol(z)
  z <- matrix(ncol = p, nrow = N)
  Omega_t <- array(dim = c(p,p,N))
  for (t in 1:N) {
    Omega_t[ , , t] <- (diag(p) + t * Epsilon) %*% Omega
    z[t, ] <- x[t, ] %*% solve(t(Omega_t[, , t]))
  }
  return(z)
}

testApproxJD <- function(){
  acf(x, lag.max = 12, type = "covariance")$acf -> covxxx
  covs <- array(dim = c(3,3, 13))
  for(i in 1:13) covs[,,i] <- covxxx[i,,]
  print("Shoud be permutation matrix")
  approxJD(covs)$W %*% solve(SOBI(x)$W) %>% round(3) %>% print
}

nearestSPD <- function(X){
  if(length(dim(X)) != 2 ){
    stop("only works for matrix")
  } else if(dim(X)[1] != dim(X)[2]){
    stop("only works for square matrix")
  }
  
  RES <- Matrix::nearPD(X, doSym = TRUE)
  
  list(mat = RES$mat, normF = RES$normF)
}


# algorithm ---------------------------------------------------------------


cov_sep_vec <- function(x, lags = 6, quadratic = TRUE, fix_symmetry = TRUE, verbose = TRUE){
  T <- nrow(x); p <- ncol(x)
  
  # lags accept two types | we need 0 for whitening
  if(length(lags) == 1) lags <- 0:lags 
  if(!(0 %in% lags)) lags <- c(0, lags)
  
  # each lags
  Beta_1 <- Beta_2 <- Beta_3 <- array(dim = c(p, p, length(lags) ))
  H_vec <- S_vec <- lm_res <- list()
  for(i in 1:length(lags)) {
    l <- lags[i]
    
    seq <- 1:(T-l)
    if(quadratic)  H <- matrix(c(rep(1, T-l), seq, seq * (seq + l)), ncol = 3)
    if(!quadratic) H <- matrix(c(rep(1, T-l), seq), ncol = 2)
    H_vec[[i]] <- H %x% diag(rep(1, p^2))
    
    S_t <- lapply(1:(T-l), function(tt) matrix(x[tt, ], ncol = 1) %*% matrix(x[tt + l, ], nrow = 1))
    S_vec[[i]] <- unlist(lapply(S_t, as.vector))
    
    lm_res[[i]] <- lm(S_vec[[i]] ~ H_vec[[i]] - 1)
    
    Beta_col_vecs  <- matrix(lm_res[[i]]$coefficients,  ncol = ifelse(quadratic, 3, 2))
    Beta_1[ , , i] <- matrix(Beta_col_vecs[ ,1], ncol = p)
    Beta_2[ , , i] <- matrix(Beta_col_vecs[ ,2], ncol = p)
    if(quadratic) Beta_3[ , , i] <- matrix(Beta_col_vecs[ ,3], ncol = p)
  }
  
  if(fix_symmetry){
    for(i in 1:length(lags)) {
      Beta_2[ , , i] <- 0.5 * (Beta_2[ , , i] + t(Beta_2[ , , i])) # apply symmetry fix
      Beta_3[ , , i] <- 0.5 * (Beta_3[ , , i] + t(Beta_3[ , , i])) # apply symmetry fix
    }
  }
  
  if(verbose){
    cat("R_squared in Cov-Separation")
    print(unlist(lapply(lm_res, function(res) summary(res)$r.squared)))
  }
  
  list(beta_1 = Beta_1, beta_2 = Beta_2, beta_3 = Beta_3, lags = lags)
}

approxJD <- function(covs, method = "frjd"){
  # ApproxJD, Miettinen, Nordhausen & Taskinen (2017)
  # param covs as in JADE | in the form of dim = p,p,lag.max+1 | lag 0 at index 1

  require(JADE)
  
  if(dim(covs)[1] != dim(covs)[2]) 
    stop('Autocovariance matrices should be in dim(p,p,k+1)')

    # here, for sure 0 should be removed
  lag_len <- dim(covs)[3] - 1 
  
    # whitening using var^{-1/2}, through eigen calculation
  nearestDist <- 0
  EVD <- eigen(covs[ , , 1], symmetric = TRUE)
  if (min(EVD$values) < 0) {
    # not semi-definite?
    spd <- nearestSPD(covs[ , , 1])
    EVD <- eigen(spd$mat, symmetric = TRUE)
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

solve_main <- function(cov_sep_res){
  require(matrixcalc)
  beta_1 <- cov_sep_res$beta_1
  beta_2 <- cov_sep_res$beta_2
  lags   <- cov_sep_res$lags
  p      <- dim(beta_1)[1]
  cov_hat<- array(dim = dim(beta_1))
    # estimates for Omega * Lambda_l * Omega
  for(i in 1:length(lags)) cov_hat[,,i] <- 0.5 * (beta_1[,,i] + t(beta_1[,,i]) - lags[i] * beta_2[,,i])
  
    # Joint Diagnolization
  JD_res <- approxJD(cov_hat)
  Omega_hat  <- solve(JD_res$W)
  Lambda_hat <- JD_res$D
    
    # Find Epsilon
  I <- diag(rep(1,p))
  H <- beta_2_vec <- NULL
  for(j in 1:dim(Lambda_hat)[3]){
    OLO_hat_j    <- Omega_hat %*% Lambda_hat[,,j] %*% t(Omega_hat)
    H_j          <- OLO_hat_j %x% I + (I %x% OLO_hat_j) %*% K.matrix(p)
    beta_2_vec_j <- matrix(c(beta_2[,,j]), nrow = p^2)
    
    H            <- rbind(H, H_j)
    beta_2_vec   <- rbind(beta_2_vec, beta_2_vec_j)
  } 
  
  lm3_res <- lm(beta_2_vec ~ H - 1)
  Epsilon_hat <- matrix(lm3_res$coefficients, nrow = p)
  
    ### info ### 
  cat("R_squared in Solving Omega: ", summary(lm3_res)$r.squared, "\n")
  
  list(W = solve(Omega_hat), Omega_hat = Omega_hat, Epsilon_hat = Epsilon_hat)
}

solve_alt <- function(cov_sep_res){
  # HUOM: beta_2[p, p, lags + 1] ; D[p, p, lags] due to covariance
  require(matrixcalc)
  beta_2 <- cov_sep_res$beta_2
  beta_3 <- cov_sep_res$beta_3
  p      <- dim(beta_2)[1]
  
  # Joint Diagnolization
  JD_res <- approxJD(beta_3)
  
  # get Epsion and Omega
  EpsilonOmega_hat <- solve(JD_res$W)
  D <- JD_res$D
  I <- diag(rep(1, p)) # identity matrix at proper dim
  H <- beta_2_vec <- NULL       # H and Y are stacked matrix
  for(i in 1:dim(D)[3]){
    A            <- EpsilonOmega_hat %*% D[,,i]
    H_i          <- (I %x% A) %*% K.matrix(p) + (A %x% I)
    beta_2_vec_i <-  matrix(c(beta_2[,,i]), nrow = p^2)
    
    H            <- rbind(H, H_i)
    beta_2_vec   <- rbind(beta_2_vec, beta_2_vec_i)
  }
  
  lm2_res <- lm(beta_2_vec ~ H - 1)
  
  ### info ###
  cat("R_squared in Solving Omega: ", summary(lm2_res)$r.squared, "\n")
  
  Omega_hat <- matrix(lm2_res$coefficients, nrow = p, byrow = FALSE)
  
  Epsilon_hat <- EpsilonOmega_hat %*% Omega_hat
  
  list(W = solve(Omega_hat), Omega_hat = Omega_hat, Epsilon_hat = Epsilon_hat)
}


# packed functions --------------------------------------------------------

tvsobi011 <- function(x, lags = 12) solve_main(cov_sep_vec(x, lags, TRUE,  TRUE))
tvsobi010 <- function(x, lags = 12) solve_main(cov_sep_vec(x, lags, TRUE,  FALSE))
tvsobi001 <- function(x, lags = 12) solve_main(cov_sep_vec(x, lags, FALSE, TRUE))
tvsobi000 <- function(x, lags = 12) solve_main(cov_sep_vec(x, lags, FALSE, FALSE))
tvsobi2   <- function(x, lags = 12) solve_alt(cov_sep_vec(x, lags))
