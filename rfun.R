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
  N <- nrow(x);  p <- ncol(x)
  z <- matrix(ncol = p, nrow = N)
  Omega_t <- array(dim = c(p,p,N))
  for (t in 1:N) {
    Omega_t[ , , t] <- (diag(p) + t * Epsilon_hat) %*% Omega_hat
    z[t, ] <- x[t, ] %*% solve(t(Omega_t[, , t]))
  }
  return(as.ts(z))
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

cov_sep <- function(x, lags = 6, quadratic = TRUE, fix_symmetry = TRUE, verbose = FALSE){
    
  N <- nrow(x); p <- ncol(x)
  
  # lags accept two types | we need 0 for whitening
  if(length(lags) == 1) lags <- 0:lags 
  if(!(0 %in% lags)) lags <- c(0, lags)
  
  # looping for each lags and element-wise matrix
  Beta_1 <- Beta_2 <- Beta_3 <- r_squared <- array(dim = c(p,p, length(lags)))
  
  for(id in 1:length(lags)) { for (i in 1:p) { for (j in 1:p) {
    
    l <- lags[id]
    
    # build y
    y <- numeric(N-l)
    for (t in 1:(N-l)){
      cov_x_t_l <- matrix(x[t, ], nrow = p) %*% matrix(x[t + l, ], ncol = p)
      y[t] <- cov_x_t_l[i, j]
    }
    # build H, except for 1 columns 
    H2 <-  1:(N-l)
    H3 <- (1:(N-l))^2 + (1:(N-l)) * l    # "Correct"
    lm1_data = data.frame(y, H2, H3)
    
    # regression for betas in alternative ways
    if( quadratic) lm1_res <- lm(y ~ H2 + H3, lm1_data)
    if(!quadratic) lm1_res <- lm(y ~ H2, lm1_data)
    
    # get result
    r_squared[i,j,id] <- summary(lm1_res)$r.squared
    Beta_1[i,j,id]    <- lm1_res$coefficients[1]
    Beta_2[i,j,id]    <- lm1_res$coefficients[2]
    Beta_3[i,j,id]    <- lm1_res$coefficients[3]
  }}} # end loop for i,j,l
  
  # fix symmetry for beta_2 and beta_3
  if(fix_symmetry){
    for(i in 1:length(lags)) {
      Beta_2[ , , i] <- 0.5 * (Beta_2[ , , i] + t(Beta_2[ , , i])) # apply symmetry fix
      Beta_3[ , , i] <- 0.5 * (Beta_3[ , , i] + t(Beta_3[ , , i])) # apply symmetry fix
    }
  }
  
  if(verbose){
    print(verbose)
    cat("R_squared in Cov-Separation")
    lapply(r_squared, mean) %>% unlist %>% print
  }
  
  method <- paste(ifelse(quadratic, "Quadratic", "Linear"),
                  ifelse(fix_symmetry, "Symmetric", "Non-Symmetric"))
  list(beta_1 = Beta_1, beta_2 = Beta_2, beta_3 = Beta_3, lags = lags, x = x, method = method)
}

cov_sep_vec <- function(x, lags = 6, quadratic = TRUE, fix_symmetry = TRUE, verbose = FALSE){
  N <- nrow(x); p <- ncol(x)
  
  # lags accept two types | we need 0 for whitening
  if(length(lags) == 1) lags <- 0:lags 
  if(!(0 %in% lags)) lags <- c(0, lags)
  
  # each lags
  Beta_1 <- Beta_2 <- Beta_3 <- array(dim = c(p, p, length(lags) ))
  H_vec <- S_vec <- lm_res <- list()
  for(i in 1:length(lags)) {
    l <- lags[i]
    
    seq <- 1:(N-l)
    if(quadratic)  H <- matrix(c(rep(1, N-l), seq, seq * (seq + l)), ncol = 3)
    if(!quadratic) H <- matrix(c(rep(1, N-l), seq), ncol = 2)
    H_vec[[i]] <- H %x% diag(rep(1, p^2))
    
    S_t <- lapply(1:(N-l), function(tt) matrix(x[tt, ], ncol = 1) %*% matrix(x[tt + l, ], nrow = 1))
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
    print(verbose)
    cat("R_squared in Cov-Separation")
    print(unlist(lapply(lm_res, function(res) summary(res)$r.squared)))
  }
  
  method <- paste(ifelse(quadratic, "Quadratic", "Linear"),
                  ifelse(fix_symmetry, "Symmetric", "Non-Symmetric"))
  list(beta_1 = Beta_1, beta_2 = Beta_2, beta_3 = Beta_3, lags = lags, x = x, method = method)
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
    # warning(paste("covariance is NOT semi-definite, applying Nearest SPD;",
    #               "Frobenius norm:", spd$normF))
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
  method <- cov_sep_res$method
  p      <- dim(beta_1)[1]
  cov_hat<- array(dim = dim(beta_1))
    # estimates for Omega * Lambda_l * Omega
  for(i in 1:length(lags)) cov_hat[,,i] <- 0.5 * (beta_1[,,i] + t(beta_1[,,i]) - lags[i] * beta_2[,,i])
  
    # Joint Diagnolization
  JD_res <- approxJD(cov_hat)
  Omega_hat  <- solve(JD_res$W)
  Lambda_hat <- JD_res$D
  if(!(JD_res$nearestDist == 0)) method <- paste(method, "NearestSPD")
    
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
  
  ### info ### cat("R_squared in Solving Omega: ", summary(lm3_res)$r.squared, "\n")
  
  S_hat <- tvunmix(cov_sep_res$x, Omega_hat, Epsilon_hat)
  res <- list(W = solve(Omega_hat), k = lags[-1], method = paste("LTV-SOBI", method), S = S_hat, Omega_hat = Omega_hat, Epsilon_hat = Epsilon_hat)
  attr(res, "class") <- "tvbss"
  return(res)
}

solve_alt <- function(cov_sep_res){
  # HUOM: beta_2[p, p, lags + 1] ; D[p, p, lags] due to covariance
  require(matrixcalc)
  beta_2 <- cov_sep_res$beta_2
  beta_3 <- cov_sep_res$beta_3
  method <- cov_sep_res$method
  lags   <- cov_sep_res$lags
  p      <- dim(beta_2)[1]
  
  # Joint Diagnolization
  JD_res <- approxJD(beta_3)
  if(!(JD_res$nearestDist == 0)) method <- paste(method, "NearestSPD")
  
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
  
  ### info ###  cat("R_squared in Solving Omega: ", summary(lm2_res)$r.squared, "\n")
  
  Omega_hat <- matrix(lm2_res$coefficients, nrow = p, byrow = FALSE)
  
  Epsilon_hat <- EpsilonOmega_hat %*% Omega_hat
  
  S_hat <- tvunmix(cov_sep_res$x, Omega_hat, Epsilon_hat)
  res <- list(W = solve(Omega_hat), k = lags[-1], method = paste("LTV-SOBI-alt", method), S = S_hat, Omega_hat = Omega_hat, Epsilon_hat = Epsilon_hat)
  attr(res, "class") <- "tvbss"
  return(res)
}

# performance measurement -------------------------------------------------

# JADE::SIR: get correlation; find correct permutation (max at row/col); get mean;
# rstudioapi::viewer("https://hal.inria.fr/inria-00544230/document")

get_WA_for_SIR <- function(res, N, Omega, Epsilon){
  Omega_hat <- res$Omega_hat
  Epsilon_hat <- res$Epsilon_hat

  p <- ncol(Omega_hat); 
  WA <- WA_r <- WA_r2 <- array(dim = c(p,p,N))
  for (t in 1:N) {
    Omega_hat_t <- (diag(p) + t * Epsilon_hat) %*% Omega_hat
    Omega_t     <- (diag(p) + t * Epsilon) %*% Omega
    WA[ , , t]  <- solve(Omega_hat_t) %*% Omega_t
    
    # find the best permutation here
    for (i in 1:p) {
      cur_col <- abs(WA[i, , t])
      id      <- which(cur_col == max(cur_col))[1]
      WA_r[i,  i, t] <- cur_col[ id]
      WA_r[i, -i, t] <- cur_col[-id]
    }
    WA_r2[ , , t] <- (WA_r[ , , t])^2 # Hadamard element wise 
  }
  
  list(WA = WA, WA_rotated = WA_r, WA_rotated_sq = WA_r2)
}

get_SIR_from_WA <- function(WA_rotated){
  
  N <- dim(WA_rotated)[3]
  p <- dim(WA_rotated)[1]
  ratio <- sum_diag <- sum_off <- numeric(N)
  
  for (i in 1:N){
    items_diag  <- diag(WA_rotated[,,i])
    items_off   <- as.vector(WA_rotated[,,i] - diag(items_diag, p)) # diag 0 is ok
    sum_diag[i] <- sum(items_diag)
    sum_off [i] <- sum(items_off )
    ratio[i]    <- sum_diag[i] / sum_off[i]
  }
  
  ratio_vec <- 10 * log10(ratio)
  ratio_new <- 10 * log10( sum(sum_diag) / sum(sum_off))
  # list(ratio1 = mean(ratio_vec), ratio2 = ratio_new)
  
  ratio_new
}

SIR_all <- function(bss_res, Omega, Epsilon, S){
  # input true Omega Epsilon and Source
  N <- nrow(bss_res$S)
  p <- ncol(bss_res$W)
  Matrix_0 <- matrix(rep(0, p^2), ncol = p)
  
  if(class(bss_res) == "tvbss"){
    WA  <- get_WA_for_SIR(bss_res, N, Omega, Epsilon)
    Omega_hat   <- bss_res$Omega_hat
    Epsilon_hat <- bss_res$Epsilon_hat
  }
  else if(class(bss_res) == "bss") {
    WA  <- get_WA_for_SIR(list(Omega_hat = solve(bss_res$W), Epsilon_hat = Matrix_0), N, Omega, Epsilon)
    Omega_hat   <- solve(bss_res$W)
    Epsilon_hat <- Matrix_0
  }
  else stop("incompatible input")
  
  MD_all_t <- numeric(N)
  for (t in 1:N) {
    MD_all_t[t] <- MD((diag(p) + t * Epsilon) %*% Omega, (diag(p) + t * Epsilon_hat) %*% Omega_hat)
  }
  
  
  list(method = bss_res$method,
       SIR_diag_sq  = get_SIR_from_WA(WA$WA_rotated_sq),
       SIR_diag_abs = get_SIR_from_WA(WA$WA_rotated),
       SIR_original = SIR(S, bss_res$S),
       MD_0         = MD(Omega_hat, Omega),
       MD_mean      = mean(MD_all_t),
       N=N, p=p)
}



# packed functions --------------------------------------------------------

ltvsobi2   <- function(x, use_vec = FALSE, lags = 12) {
  x <- scale(x, center = TRUE, scale = FALSE)
  if(!use_vec) return(solve_alt(cov_sep(x, lags)))
  if(use_vec)  return(solve_alt(cov_sep_vec(x, lags)))
}

ltvsobi <- function(x, lags = 12, quadratic = TRUE, fix_symmetry = TRUE, use_vec = TRUE, verbose = FALSE) {
  x <- scale(x, center = TRUE, scale = FALSE)
  if( use_vec) return(solve_main(cov_sep_vec(x, lags, quadratic = quadratic, fix_symmetry = fix_symmetry, verbose = verbose)) )
  if(!use_vec) return(solve_main(cov_sep    (x, lags, quadratic = quadratic, fix_symmetry = fix_symmetry, verbose = verbose)) )
}
