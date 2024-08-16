# initalizations
source('00_system.R')

args <- commandArgs(trailingOnly = TRUE)
loadstring <- as.character(args[1])
print(loadstring)

tt <- sum(.Internal(gc(FALSE, TRUE, TRUE))[13:14]) # memory report

elapsed_time <- system.time({
  
  load(file = paste0(CommonPath, '/fields/', loadstring))
  
  if(this_field@gaussian){
    
    # Estimate betas
    distmat <- as.matrix(dist(this_field@locs))
    Sigma_oracle <- as.matrix(cov.mat2(h = distmat, theta = this_field@truetheta))
    z <- this_field@fulldata
    X <- cbind(rep(1, length(this_field@fulldata)), cos(this_field@locs))
    n <- length(z)
    colnames(X) <- NULL
    
    xsigma <- t(solve(Sigma_oracle, X))
    right_side <- xsigma %*% as.matrix(z, ncol = 1)
    beta0_init <- solve(xsigma %*% X, right_side)
    
    z_nomean <- z - ( X %*% beta0_init )
    
    theta_lower <- c(1e-8, 1e-8, 1e-8)
    theta_upper <- c(3 * sqrt(2), 20, 1)
    theta_0 <-  c(0.1, 0.9*var(z_nomean), 0.1*var(z_nomean))
    
    output <- GeoFit(data = z_nomean, coordx = this_field@locs, corrmodel = 'Matern', 
                     model = 'Gaussian', maxdist = ifelse(this_field@gridded, 0.04, 0.05 ), likelihood = 'Conditional',  type = 'Pairwise',
                     optimizer = 'L-BFGS-B', parallel = TRUE,
                     
                     lower = list(scale = theta_lower[1], sill = theta_lower[2], nugget = theta_lower[3]), 
                     
                     upper = list(scale = theta_upper[1], sill = theta_upper[2], nugget = theta_upper[3]),
                     
                     start = list(scale = theta_0[1], sill = theta_0[2], nugget = theta_0[3]), 
                     
                     fixed = list(smooth = this_field@truetheta[3], mean = 0))
    
  } else {
    
    X <- cbind(rep(1, length(this_field@fulldata)), cos(this_field@locs))
    
    theta_lower <- c(1e-8, 1e-8, 1e-8, -5, 1e-8, -20, - 20, -20)
    theta_upper <- c(3 * sqrt(2), 20, 1, 5, 5, 20, 20, 20)
    theta_0 <- c(0.1, 0.9 * var(this_field@fulldata), 0.1 * var(this_field@fulldata), 0, 1, mean(this_field@fulldata), 0, 0)
    
    output <- GeoFit(data = this_field@fulldata, coordx = this_field@locs, corrmodel = 'Matern', 
                     model = 'SinhAsinh', maxdist = ifelse(this_field@gridded, 0.04, 0.05 ), likelihood = 'Conditional',  type = 'Pairwise',
                     optimizer = 'L-BFGS-B', parallel = TRUE, X = X,
                     
                     lower = list(scale = theta_lower[1], sill = theta_lower[2],
                                  nugget = theta_lower[3], skew = theta_lower[4], 
                                  tail = theta_lower[5], mean = theta_lower[6], 
                                  mean1 = theta_lower[7], mean2 = theta_lower[8]), 
                     
                     upper = list(scale = theta_upper[1], sill = theta_upper[2],
                                  nugget = theta_upper[3], skew = theta_upper[4], 
                                  tail = theta_upper[5], mean = theta_upper[6], 
                                  mean1 = theta_upper[7], mean2 = theta_upper[8]),
                     
                     start = list(scale = theta_0[1], sill = theta_0[2],
                                  nugget = theta_0[3], skew = theta_0[4], 
                                  tail = theta_0[5], mean = theta_0[6], 
                                  mean1 = theta_0[7], mean2 = theta_0[8]), 
                     
                     fixed = list(smooth = this_field@truetheta[3])
                     
                     )
    
  }
  
})[3]

ifelse(this_field@gaussian, betas_hat <- beta0_init,
       betas_hat <- output$param[grep('mean', names(output$param))])

memory_used <- sum(.Internal(gc(FALSE, FALSE, TRUE))[13:14]) - tt

this_estimates <- new('composite_estimate',
                      specsnr = this_field@specsnr,
                      fielditer = this_field@fielditer,
                      time = elapsed_time,
                      optimout = output,
                      betas_hat = betas_hat,
                      memory_used = memory_used,
                      loadstring = loadstring)

save(this_estimates, file = paste0(CommonPath, '/estimates/estimates.composite.', this_field@specsnr, '.iter', this_field@fielditer, '.RData'))