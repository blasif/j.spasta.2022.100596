# initalizations
source('00_system.R')

args <- commandArgs(trailingOnly = TRUE)
loadstring <- as.character(args[1])
print(loadstring)

# estimation

tt <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14]) # memory report

elapsed_time <- system.time({
  
  # get the data and simulationspecs
  load(file = paste0(CommonPath, '/fields/', loadstring))
  delta <- get_range_help(this_field@truetheta)/2
  memory_values <- get_memory_values(this_field)
  
  # Covariance function taperfunction chosen only on true smoothness
  
  Covariance <- function( h, theta, delta, ... ) {
    if (this_field@truetheta[3] <= 0.5){
      return(cov.mat2( h, theta, ... ) * cov.sph( h, c(delta, 1), ...))
    } else{
      if (this_field@truetheta[3] <= 1.5){
        return(cov.mat2( h, theta, ... ) * cov.wend1( h, c(delta, 1), ...))
      } else return(cov.mat2( h, theta, ... ) * cov.wend2( h, c(delta, 1), ...))
    }  
  }
  
  options(spam.nearestdistnnz = c(memory_values[1]/2 , 1)) # optimized
  
  if(this_field@gaussian){
    
    # Estimate betas
    distmat <- nearest.dist(this_field@locs, delta = delta, upper = NULL)
    Sigma_oracle <- Covariance(h = distmat, theta = this_field@truetheta, delta = delta)
    Rstruct <- chol(Sigma_oracle, memory = list(nnzR = memory_values[2] / 2)) # optimized
    z <- this_field@fulldata
    X <- cbind(rep(1, length(this_field@fulldata)), cos(this_field@locs))
    n <- length(z)
    colnames(X) <- NULL
    xsigma <- t(solve(Sigma_oracle, X))
    right_side <- xsigma %*% as.matrix(z, ncol = 1)
    beta0_init <- solve(xsigma %*% X, right_side)
    
    z_nomean <- z - ( X %*% beta0_init )
    
    neg2loglikelihood <- function(theta) {
      Sigma <- Covariance(h = distmat,
                          theta = c(theta[1:2], this_field@truetheta[3], theta[3]),
                          delta = delta)
      cholS <- chol(Sigma, Rstruct = Rstruct)
      logdet <- sum(log(diag(cholS)))
      
      return(n * log(2 * pi) + 2 * logdet + sum(z_nomean * backsolve(cholS,
                                                                     forwardsolve(cholS, z_nomean, transpose = TRUE, upper.tri = TRUE),
                                                                     n)))
    }
    
    theta_lower <- c(1e-8, 1e-8, 1e-8)
    theta_upper <- c(Inf, 20, 1)
    theta_0 <-  c(0.1, 0.9*var(z_nomean), 0.1*var(z_nomean))
    
    cl <- makeCluster(ifelse(ncores > length(theta_0) * 2 + 1, length(theta_0) * 2 + 1, 1 + floor(ncores/3)),
                      useXDR = F)
    setDefaultCluster(cl = cl)
    clusterExport(cl, c('lib.loc'))
    clusterEvalQ(cl, .libPaths(new = lib.loc))
    clusterEvalQ(cl, library('spam', lib.loc = lib.loc))
    clusterEvalQ(cl, library('spam64', lib.loc = lib.loc))
    clusterExport(cl, c('neg2loglikelihood', 'z_nomean', 'distmat', 'Rstruct', 'n', 'delta', 'Covariance', 
                        'this_field', 'cov.mat2'))
    
    output <- optimParallel::optimParallel(fn = neg2loglikelihood,
                                           method = 'L-BFGS-B',
                                           lower = theta_lower,
                                           upper = theta_upper,
                                           par = theta_0,
                                           parallel = list(forward = FALSE,
                                                           loginfo = TRUE),
                                           control = list(factr = 1e8, # approx. 1e-8
                                                          trace = TRUE, maxit = 1000))
    stopCluster(cl)
    
  } else {
    
    distmat <- nearest.dist(this_field@locs, delta = delta, upper = NULL)
    Sigma_oracle <- Covariance(h = distmat, theta = this_field@truetheta, delta = delta)
    Rstruct <- chol(Sigma_oracle, memory = list(nnzR = memory_values[2] / 2))
    z <- this_field@fulldata
    X <- cbind(rep(1, length(this_field@fulldata)), cos(this_field@locs))
    n <- length(z)
    colnames(X) <- NULL
    
    neg2logliksasbetas <- function(theta){
      
      Sigma <- Covariance(h = distmat,
                          theta = c(theta[1:2], this_field@truetheta[3], theta[3]),
                          delta = delta)
      
      corr_matrixd <- Sigma/Sigma[1, 1]
      cholS <- chol(corr_matrixd, Rstruct = Rstruct)
      logdet <- 2 * sum(log(diag(cholS)))
      stdata <- (z - X %*% matrix(theta[6:8], nrow = 3)) / sqrt(theta[2])
      skew <- theta[4]
      kurt <- theta[5]
      
      Z <- sinh(kurt * asinh(stdata) - skew)
      C <- sqrt(1+Z^2)
      
      return( - 2 * (n * log(kurt) - (n/2) * log(2 * pi * (theta[2])) - 0.5 * logdet + 
                       sum( log(C) - 0.5 * log(1 + stdata^2))  - 0.5 * sum(Z * backsolve(cholS,
                                                                                         forwardsolve(cholS, Z, transpose = TRUE, upper.tri = TRUE),
                                                                                         n)) )
      )
      
    }
    
    theta_lower <- c(1e-8, 1e-8, 1e-8, -5, 1e-8, -20, - 20, -20)
    theta_upper <- c(Inf, 20, 1, 5, 5, 20, 20, 20)
    theta_0 <- c(0.1, 0.9 * var(this_field@fulldata), 0.1 * var(this_field@fulldata), 0, 1, mean(this_field@fulldata), 0, 0)
    
    cl <- makeCluster(ifelse(ncores > length(theta_0)*2 +1, length(theta_0)*2 + 1, 1 + floor(ncores/3)),
                      useXDR = F)
    setDefaultCluster(cl = cl)
    clusterExport(cl, c('lib.loc'))
    clusterEvalQ(cl, .libPaths(new = lib.loc))
    clusterEvalQ(cl, library('spam', lib.loc = lib.loc))
    clusterEvalQ(cl, library('spam64', lib.loc = lib.loc))
    clusterExport(cl, c('neg2logliksasbetas', 'z', 'distmat', 'Rstruct', 'n', 'X', 'this_field', 'cov.mat2', 'Covariance', 'delta'))
    
    output <- optimParallel::optimParallel(fn = neg2logliksasbetas,
                                           method = 'L-BFGS-B',
                                           lower = theta_lower,
                                           upper = theta_upper,
                                           par = theta_0,
                                           parallel = list(forward = FALSE,
                                                           loginfo = TRUE),
                                           control = list(factr = 1e8, # approx. 1e-8
                                                          trace = TRUE, maxit = 1000))
    stopCluster(cl)
    
  }
  
})[3]

ifelse(this_field@gaussian, betas_hat <- beta0_init,
       betas_hat <- output$par[6:8])

memory_used <- sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt

this_estimates <- new('taper_estimate',
                      specsnr = this_field@specsnr,
                      fielditer = this_field@fielditer,
                      delta = delta,
                      time = elapsed_time,
                      optimout = output,
                      betas_hat = betas_hat,
                      memory_used = memory_used,
                      loadstring = loadstring)

save(this_estimates, file = paste0(CommonPath, '/estimates/estimates.taper.', this_field@specsnr, '.iter', this_field@fielditer, '.RData'))