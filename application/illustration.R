# illustration

rm(list = ls())
source('00_system.R')
load(paste0(CommonPath, '/application/ClimChange.RData'))

# load data and spam options
{
  delta <- (tas2 - tas1)[, c(2:191)]
  lati <- lati[c(2:191)]
  image.plot(long, lati, delta)
  map(add = TRUE, wrap = c(0, 360))
  lonlat <- expand.grid(lon = long-180, lati = lati)
  options(spam.nearestdistnnz = c(15326784 * 1.10, 1)) # optimized
}

# compute mean part
{
  y <- c(delta)
  xlon <- pi*(lonlat[, 1]-180)/180
  xlat <- pi*lonlat[, 2]/180
  model_test <- gam(y ~ s(sin(xlat), bs = 'cr')+ s(cos(xlat), bs = 'cr'))
  deltamean <- matrix( fitted(model_test), nrow = 288, ncol = 190)
  image.plot(long, lati, deltamean)
  deltaresid <- matrix( resid(model_test), nrow = 288, ncol = 190)
  image.plot(long, lati, deltaresid)
}

# Cov. parameters estimation ----

delta_taper <- 6
# compute tapered great circle distances

# Tapering -----

Covariance <- function( h, theta, delta_taper, ... ) {
  return(cov.wend1_1(h, theta, ...) * cov.wend2_1(h, c(delta_taper, 1, 0), ...))
}

distmat <- nearest.dist(lonlat, method = 'greatcircle', upper = NULL, delta = delta_taper)
covmiss <- Covariance(distmat, theta = c(5, 1, 0), delta_taper = delta_taper)
Rstruct <- chol(covmiss)

y <- as.vector(deltaresid)
n <- length(y)

neg2logliksas0_taper <- function(theta){
  
  Sigma <- Covariance(h = distmat,
                      theta = c(theta[1:3]),
                      delta_taper = delta_taper)
  
  corr_matrixd <- Sigma/Sigma[1, 1]
  
  cholS <- chol(corr_matrixd, Rstruct = Rstruct)
  
  logdet <- 2 * sum(log(diag(cholS)))
  
  stdata <- y / sqrt(theta[2])
  
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

cl <- makeCluster(ifelse(ncores > 5*2+1,5*2+1, 1 + sqrt(ncores)/3),
                  useXDR = F)
setDefaultCluster(cl = cl)
clusterExport(cl, c('lib.loc'))
clusterEvalQ(cl, .libPaths(new = lib.loc))
clusterEvalQ(cl, library('spam', lib.loc = lib.loc))
clusterEvalQ(cl, library('spam64', lib.loc = lib.loc))
clusterExport(cl, c('neg2logliksas0_taper', 'y', 'distmat', 'Rstruct', 'n', 'delta_taper', 'Covariance', 
                    'cov.wend1_1', 'cov.wend2_1'))

theta_lower <- c(1e-8, 1e-8, 1e-8, -5, 1e-8)
theta_upper <- c(50, 5*var(y), 0.5, 5, 5)
theta_0 <- c(25, 0.9 * var(y), 0.1 * var(y), 0, 1)

output <- optimParallel::optimParallel(fn = neg2logliksas0_taper,
                                       method = 'L-BFGS-B',
                                       lower = theta_lower,
                                       upper = theta_upper,
                                       par = theta_0,
                                       parallel = list(forward = FALSE,
                                                       loginfo = FALSE),
                                       control = list(factr = 1e8, # approx. 1e-8
                                                      trace = FALSE, maxit = 1000), 
                                       hessian = F)
stopCluster(cl)

output_taper <- output
save(output_taper, file = paste0(CommonPath, '/application/Output_taper.RData'))
load(paste0(CommonPath, '/application/Output_taper.RData'))

# compute the hessian w/ 40k obs

set.seed(1005)
sample_test <- sample(1:dim(lonlat)[1], 40000)
distmat <- nearest.dist(lonlat[sample_test, ], method = 'greatcircle', upper = NULL, delta = delta_taper)
covmiss <- Covariance(h = distmat,
                      theta = output_taper$par,
                      delta_taper = delta_taper)
Rstruct <- chol(covmiss)
summary(Rstruct)
y <- as.vector(deltaresid)[sample_test]
hessian_taper <- hessian(neg2logliksas0_taper, x = output_taper$par)

sqrt(diag(solve(hessian_taper))) # 11.212679524  0.002166104  0.167379648  0.002347452  0.007525205

save(hessian_taper, file = paste0(CommonPath, '/application/hessian_taper.RData'))
load(paste0(CommonPath, '/application/hessian_taper.RData'))

# Deliberate Misspecification ------

distmat <- nearest.dist(lonlat, method = 'greatcircle', upper = NULL, delta = delta_taper)
covmiss <- cov.wend2_1(h = distmat,
                       theta = c(delta_taper, 1, 0))
Rstruct <- chol(covmiss)

y <- as.vector(deltaresid)
n <- length(y)

neg2logliksas0_dm <- function(theta){
  
  Sigma <- cov.wend2_1(h = distmat,
                       theta = c(delta_taper, theta[1:2]))
  
  corr_matrixd <- Sigma/Sigma[1, 1] 
  
  cholS <- chol(corr_matrixd, Rstruct = Rstruct)
  
  logdet <- 2 * sum(log(diag(cholS)))
  
  stdata <- y / sqrt(theta[1])
  
  skew <- theta[3]
  kurt <- theta[4]
  
  Z <- sinh(kurt * asinh(stdata) - skew)
  C <- sqrt(1+Z^2)
  
  return( - 2 * (n * log(kurt) - (n/2) * log(2 * pi * (theta[1])) - 0.5 * logdet + 
                   sum( log(C) - 0.5 * log(1 + stdata^2))  - 0.5 * sum(Z * backsolve(cholS,
                                                                                     forwardsolve(cholS, Z, transpose = TRUE, upper.tri = TRUE),
                                                                                     n)) )
  )
  
}

cl <- makeCluster(ifelse(ncores > 4*2+1,4*2+1, 1 + sqrt(ncores)/3),
                  useXDR = F)
setDefaultCluster(cl = cl)
clusterExport(cl, c('lib.loc'))
clusterEvalQ(cl, .libPaths(new = lib.loc))
clusterEvalQ(cl, library('spam', lib.loc = lib.loc))
clusterEvalQ(cl, library('spam64', lib.loc = lib.loc))
clusterExport(cl, c('neg2logliksas0_dm', 'y', 'distmat', 'Rstruct', 'n', 'cov.wend2_1', 'delta_taper'))

theta_lower <- c(1e-8, 1e-8, -5, 1e-8)
theta_upper <- c(5*var(y), 0.5, 5, 5)
theta_0 <- c(0.9 * var(y), 0.1 * var(y), 0, 1)

output <- optimParallel::optimParallel(fn = neg2logliksas0_dm,
                                       method = 'L-BFGS-B',
                                       lower = theta_lower,
                                       upper = theta_upper,
                                       par = theta_0,
                                       parallel = list(forward = FALSE,
                                                       loginfo = TRUE),
                                       control = list(factr = 1e8, # approx. 1e-8
                                                      trace = TRUE, maxit = 1000), 
                                       hessian = T)
stopCluster(cl)

output_directmiss <- output

save(output_directmiss, file = paste0(CommonPath, '/application/Output_directmiss.RData'))
load(paste0(CommonPath, '/application/Output_directmiss.RData'))

# hessian dm ----

set.seed(seed_master_number + 121) # 81
sample_test <- sample(1:dim(lonlat)[1], 40000)
distmat <- nearest.dist(lonlat[sample_test, ], method = 'greatcircle', upper = NULL, delta = delta_taper)
covmiss <- cov.wend2_1(h = distmat,
                       theta = c(delta_taper, 1, 0))
Rstruct <- chol(covmiss)
summary(Rstruct)
y <- as.vector(deltaresid)[sample_test]
n <- length(y)

hessian_dm <- hessian(neg2logliksas0_dm, x = output_directmiss$par)

sqrt(diag(solve(hessian_dm))) # 0.003004825 0.032757661 0.002349393 0.002091773

save(hessian_dm, file = paste0(CommonPath, '/application/hessian_dm.RData'))
load(paste0(CommonPath, '/application/hessian_dm.RData'))

# Composite -----

y <- as.vector(deltaresid)
n <- length(y)

theta_lower <- c(1e-8, 1e-8, 1e-8, -5, 1e-8)
theta_upper <- c(40, 5*var(y), 0.5, 5, 5)
theta_0 <- c(25, 0.9 * var(y), 0.1 * var(y), 0, 1)

output <- GeoFit(data = y, coordx = lonlat, corrmodel = 'Wend2', distance = 'Geod',
                 model = 'SinhAsinh', maxdist = 6, likelihood = 'Conditional',  type = 'Pairwise',
                 optimizer = 'L-BFGS-B', parallel = TRUE,
                 
                 lower = list(scale = theta_lower[1], sill = theta_lower[2],
                              nugget = theta_lower[3], skew = theta_lower[4], 
                              tail = theta_lower[5]), 
                 
                 upper = list(scale = theta_upper[1], sill = theta_upper[2],
                              nugget = theta_upper[3], skew = theta_upper[4], 
                              tail = theta_upper[5]),
                 
                 start = list(scale = theta_0[1], sill = theta_0[2],
                              nugget = theta_0[3], skew = theta_0[4], 
                              tail = theta_0[5]), 
                 
                 fixed = list(power2 = 1.5 , mean = 0),
                 
)

output_geofit <- output
save(output_geofit, file = paste0(CommonPath, '/application/Output_geofit.RData'))
load(paste0(CommonPath, '/application/Output_geofit.RData'))

save(y, output_taper, output_directmiss, output_geofit, file = paste0(CommonPath, '/application/illustration.RData'))
