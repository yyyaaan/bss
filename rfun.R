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

# pseudo data -------------------------------------------------------------

run <- function(){
  N <- 1e3
  z <-cbind(var1 = arima.sim(list(ar=runif(1,-1,1)),N),
            var2 = arima.sim(list(ar=runif(1,-1,1)),N),
            var3 = arima.sim(list(ar=runif(1,-1,1)),N))
  Omega <- matrix(runif(9, 1, 10), ncol = 3)
  Epsilon <- 1e-4 * matrix(runif(9, 1, 10), ncol = 3)
  x <- tvmix(z, Omega, Epsilon)
  lags <- 6
  
  x %>% cov_sep %>% use_series(beta_3) %>% 
    approxJD() %>% use_series(W) %>% 
    MD(Omega %*% Epsilon) %>% print
  
  x %>% cov_sep_vec %>% use_series(beta_3) %>% 
    approxJD() %>% use_series(W) %>% 
    MD(Omega %*% Epsilon) %>% print
  
  ss <- tvsobi(x, lag.max = 12); MD(ss$W, Omega); MD(ss$Epsilon, solve(Epsilon))
  ss <- tvsobi2(x,lags = 12); MD(ss$Omega_hat %>% solve, Omega);  MD(ss$Epsilon_hat %>% solve, Epsilon)
  ss <- tvsobi3(x,lags = 12); MD(ss$Omega_hat %>% solve, Omega);  MD(ss$Epsilon_hat %>% solve, Epsilon)

  
}



# algorithm ---------------------------------------------------------------


cov_sep_vec <- function(x, lags = 6){
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
    H <- matrix(c(rep(1, T-l), seq, seq * (seq + l)), ncol = 3)
    H_vec[[i]] <- H %x% diag(rep(1, p^2))
    
    S_t <- lapply(1:(T-l), function(tt) matrix(x[tt, ], ncol = 1) %*% matrix(x[tt + l, ], nrow = 1))
    S_vec[[i]] <- unlist(lapply(S_t, as.vector))
    
    lm_res[[i]] <- lm(S_vec[[i]] ~ H_vec[[i]] - 1)
    
    Beta_col_vecs  <- matrix(lm_res[[i]]$coefficients,  ncol = 3)
    Beta_1[ , , i] <- matrix(Beta_col_vecs[ ,1], ncol = p)
    Beta_2x <- matrix(Beta_col_vecs[ ,2], ncol = p)
    Beta_3x <- matrix(Beta_col_vecs[ ,3], ncol = p)
    Beta_2[ , , i] <- 0.5 * (Beta_2x + t(Beta_2x)) # apply symmetry fix
    Beta_3[ , , i] <- 0.5 * (Beta_3x + t(Beta_3x)) # apply symmetry fix
  }
  
  
  ### info ###
  cat("R_squared in VEC Cov-Separation")
  print(unlist(lapply(lm_res, function(res) summary(res)$r.squared)))
  
  list(beta_1 = Beta_1, beta_2 = Beta_2, beta_3 = Beta_3)
}

cov_sep_stack <- function(x, lags = 6){
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
    H <- matrix(c(rep(1, T-l), seq, seq * (seq + l)), ncol = 3)
    H_vec[[i]] <- H %x% diag(rep(1, p^2))
    
    S_t <- lapply(1:(T-l), function(tt) matrix(x[tt, ], ncol = 1) %*% matrix(x[tt + l, ], nrow = 1))
    S_vec[[i]] <- unlist(lapply(S_t, as.vector))
  }
  
  S_stack <- unlist(S_vec)
  H_stack <- Matrix::bdiag(H_vec)
  H <- as.matrix(H_stack)
  res <- lm(S_stack ~ H - 1)
  summary(res)
  
  ### info ###
  cat("R_squared in VEC Cov-Separation")
  print(unlist(lapply(lm_res, function(res) summary(res)$r.squared)))
  
  list(beta_1 = Beta_1, beta_2 = Beta_2, beta_3 = Beta_3)
}

cov_sep <- function(x, lags = 6, choice = 3, fix_symmetry = TRUE){
  N <- nrow(x); p <- ncol(x)

    # lags accept two types | we need 0 for whitening
  if(length(lags) == 1) lags <- 0:lags 
  if(!(0 %in% lags)) lags <- c(0, lags)

    # looping for each lags and element-wise matrix
  empty_mat <- array(dim = c(p,p, length(lags)))
  beta_1 <- beta_2 <- beta_3 <- r_squared <- rep(list(empty_mat), 3); # array(dim = c(p,p, length(lags)))
  
  for(id in 1:length(lags)) { for (i in 1:p) { for (j in 1:p) {

    l <- lags[id]

      # build y
    y <- numeric(N-l)
    for (t in 1:(N-l)){
      cov_x_t_l <- matrix(x[t, ], nrow = 3) %*% matrix(x[t + l, ], ncol = 3)
      y[t] <- cov_x_t_l[i, j]
    }
      # build H, except for 1 columns 
    H2 <-  1:(N-l)
    H3a <- (1:(N-l))^2                    #  Yeredor
    H3b <- (1:(N-l))^2 + (1:(N-l)) * l    # "Correct"
    lm1_data = data.frame(y, H2, H3a, H3b)
    
      # regression for betas in alternative ways
    lm1_res <- list()
    lm1_res[[1]] <- lm(y ~ H2, lm1_data)
    lm1_res[[2]] <- lm(y ~ H2 + H3a, lm1_data)
    lm1_res[[3]] <- lm(y ~ H2 + H3b, lm1_data)

      # get result
    for(alt in 1:3){
      r_squared[[alt]][i,j,id] <- summary(lm1_res[[alt]])$r.squared
      beta_1[[alt]][i,j,id]    <- lm1_res[[alt]]$coefficients[1]
      beta_2[[alt]][i,j,id]    <- lm1_res[[alt]]$coefficients[2]
      beta_3[[alt]][i,j,id]    <- lm1_res[[alt]]$coefficients[3]
    }
  }}} # end loop for i,j,l
  
  # fix symmetry for beta_2 and beta_3
  if (fix_symmetry){
    for(alt in 1:3){
      for(id in 1:length(lags))
        beta_2[[alt]][,,id] = 0.5 * (beta_2[[alt]][,,id] + t(beta_2[[alt]][,,id]))
        beta_3[[alt]][,,id] = 0.5 * (beta_3[[alt]][,,id] + t(beta_3[[alt]][,,id]))
    }
  }
  
  ### info ###
  cat("R_squared in Cov-Separation: ")
  lapply(r_squared, mean) %>% unlist %>% print
  return(list(beta_1 = beta_1[[choice]], beta_2 = beta_2[[choice]], beta_3 = beta_3[[choice]]))
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

solve_omega <- function(JD_res, beta_2){
  # HUOM: beta_2[p, p, lags + 1] ; D[p, p, lags] due to covariance
  require(matrixcalc)
  EpsilonOmega_hat <- solve(JD_res$W)
  D <- JD_res$D
  p <- dim(EpsilonOmega_hat)[1]

  I <- diag(rep(1, p)) # identity matrix at proper dim
  H <- beta_2_vec <- NULL       # H and Y are stacked matrix
  for(i in 1:dim(D)[3]){
    A <- EpsilonOmega_hat %*% D[,,i]
    H_i <- (I %x% A) %*% K.matrix(p) + (A %x% I)
    beta_2_vec_i <-  matrix(c(beta_2[,,i]), nrow = p^2)
    H <- rbind(H, H_i)
    beta_2_vec <- rbind(beta_2_vec, beta_2_vec_i)
  }
  
  lm2_res <- lm(beta_2_vec ~ H - 1)
  
  ### info ###
  cat("R_squared in Solving Omega: ", summary(lm2_res)$r.squared, "\n")
  
  Omega_hat <- matrix(lm2_res$coefficients, nrow = p, byrow = FALSE)
  
  Epsilon_hat <- EpsilonOmega_hat %*% Omega_hat
  
  list(Omega_hat = Omega_hat, Epsilon_hat = Epsilon_hat)
}

tvsobi2 <- function(x, lags){
  step1 <- cov_sep(x, lags)
  step2 <- approxJD(step1$beta_3)
  solve_omega(step2, step1$beta_2)
}

tvsobi3 <- function(x, lags){
  step1 <- cov_sep_vec(x, lags)
  step2 <- approxJD(step1$beta_3)
  solve_omega(step2, step1$beta_2)
}


# imported function -------------------------------------------------------


MD_fun <- function (W.hat, A) {
  G <- W.hat %*% A
  RowNorms <- sqrt(rowSums(G^2))
  G.0 <- sweep(G, 1, RowNorms, "/")
  G.tilde <- G.0^2
  p <- nrow(A)
  Pmin <- pMatrix.min(G.tilde)
  G.tilde.p <- Pmin$A
  md <- sqrt(p - sum(diag(G.tilde.p)))/sqrt(p - 1)
  md
}

pMatrix.min <- function (A) 
{
  cost <- t(apply(A^2, 1, sum) - 2 * A + 1)
  vec <- c(solve_LSAP(cost))
  list(A = A[vec, ], pvec = vec)
}

solve_LSAP <- function (x, maximum = FALSEALSE) {
  if (!is.matrix(x) || any(x < 0)) 
    stop("x must be a matrix with nonnegative entries.")
  nr <- nrow(x)
  nc <- ncol(x)
  if (nr > nc) 
    stop("x must not have more rows than columns.")
  if (nc > nr) 
    x <- rbind(x, matrix(2 * sum(x), nc - nr, nc))
  if (maximum) 
    x <- max(x) - x
  storage.mode(x) <- "double"
  out <- .C(C_solve_LSAP, x, as.integer(nc), p = integer(nc))$p + 
    1
  out <- out[seq_len(nr)]
  class(out) <- "solve_LSAP"
  out
}


