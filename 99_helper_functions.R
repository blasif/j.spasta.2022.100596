get_range_help <- function(thetas){
  values_min <- round(0.1 * thetas[1], 3)
  values_med <- thetas[1]
  values_high <- 10 * thetas[1]
  
  to_test <- function(eff.range) {log((0.05 * thetas[2] -
                                         cov.mat(h = eff.range, theta = thetas))^2)}
  
  results <- optim(fn = to_test, par = values_med,
                   lower = values_min,
                   upper = values_high,
                   method = 'L-BFGS-B')
  
  return(results$par)
  
}

get_memory_values <- function(this_field){
  
  return(switch(as.character(this_field@nsize),
         
         '784' = switch(this_field@parsetting,
                      
                      'low' = switch(as.character(this_field@gridded),
                                     
                                     'TRUE' = c(30592, 59217),
                                     'FALSE' = c(30592, 47127)
                                     
                                     ),
                      
                      'high' = switch(as.character(this_field@gridded),
                                      
                                      'TRUE' = c(96828, 120489),
                                      'FALSE' = c(96828, 100884)
                                      
                                      
                                      )
                      
                      ),
         
         '3025' = switch(this_field@parsetting,
                       
                       'low' = switch(as.character(this_field@gridded),
                                      
                                      'TRUE' = c(462339, 1093212),
                                      'FALSE' = c(462339, 1093212)
                                      
                                      ),
                       
                       'high' = switch(as.character(this_field@gridded),
                                       
                                       'TRUE' = c(1438823, 2097923),
                                       'FALSE' = c(1500000, 2097923)
                                       
                                       )
                       
                       
                       )
         
         ) )
  
}

cov.mat2 <- function (h, theta, ..., eps = getOption('spam.eps')) 
{
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    ypos <- tmp < eps
    npos <- !ypos
    if (any(ypos)) {
      h@entries[ypos] <- theta[2]
    }
    if (any(npos)) {
      h@entries[npos] <- theta[2] * (1 - theta[4]) * (((2^(-(theta[3] - 
                                                               1)))/gamma(theta[3])) * (tmp[npos]^theta[3]) * 
                                                        besselK(tmp[npos], nu = theta[3]))
    }
  }
  else {
    tmp <- c(h)/theta[1]
    ypos <- tmp < eps
    npos <- !ypos
    if (any(ypos)) {
      h[ypos] <- theta[2]
    }
    if (any(npos)) {
      h[npos] <- theta[2] * (1 - theta[4]) * (((2^(-(theta[3] - 1)))/gamma(theta[3])) * 
                                                (tmp[npos]^theta[3]) * besselK(tmp[npos], nu = theta[3]))
    }
  }
  return(h)
}

cov.sph_1 <- function (h, theta, ..., eps = getOption('spam.eps')) 
{
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    if (theta[3] > eps) {
      ypos <- tmp <= eps
      npos <- tmp >= 1
      ypos2 <- !(ypos | npos)
      if (any(ypos)) {
        h@entries[ypos] <- theta[2]
      }
    }
    else {
      npos <- tmp >= 1
      ypos2 <- !npos
    }
    if (any(ypos2)) {
      ttmp <- tmp[ypos2]
      h@entries[ypos2] <- if (abs(theta[2] - 1) > eps) {
        theta[2] * (1 -  theta[3]) * (1 - 1.5 * ttmp + 0.5 * (ttmp * ttmp) * 
                      ttmp)
      }
      else {
        (1 - 1.5 * ttmp + 0.5 * (ttmp * ttmp) * ttmp)
      }
    }
    if (any(npos)) {
      h@entries[npos] <- 0
    }
  }
  else {
    tmp <- c(h)/theta[1]
    if (theta[3] > eps) {
      ypos <- tmp <= eps
      npos <- tmp >= 1
      ypos2 <- !(ypos | npos)
      if (any(ypos)) {
        h[ypos] <- theta[2]
      }
    }
    else {
      npos <- tmp >= 1
      ypos2 <- !npos
    }
    if (any(ypos2)) {
      ttmp <- tmp[ypos2]
      h[ypos2] <- if (abs(theta[2] - 1) > eps) 
        theta[2] * (1 - theta[3]) * (1 - 1.5 * ttmp + 0.5 * (ttmp * ttmp) * 
                      ttmp)
      else (1 - 1.5 * ttmp + 0.5 * (ttmp * ttmp) * ttmp)
    }
    if (any(npos)) {
      h[npos] <- 0
    }
  }
  return(h)
}

cov.wend1_1 <- function (h, theta, ..., eps = getOption('spam.eps')) 
{
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1
    if (any(ypos)) {
      h@entries[ypos] <- theta[2]
    }
    if (any(ypos2)) {
      h@entries[ypos2] <- theta[2] * (1 - theta[3]) * ((1 - tmp[ypos2])^4 * 
                                        (4 * tmp[ypos2] + 1))
    }
    if (any(npos)) {
      h@entries[npos] <- 0
    }
  }
  else {
    tmp <- c(h)/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1
    if (any(ypos)) {
      h[ypos] <- theta[2]
    }
    if (any(ypos2)) {
      h[ypos2] <- theta[2] * (1 - theta[3]) * ((1 - tmp[ypos2])^4 * (4 * 
                                                      tmp[ypos2] + 1))
    }
    if (any(npos)) {
      h[npos] <- 0
    }
  }
  return(h)
}

cov.wend2_1 <- function (h, theta, ..., eps = getOption('spam.eps')) 
{
  if (is.spam(h)) {
    tmp <- h@entries/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1
    if (any(ypos)) {
      h@entries[ypos] <- theta[2]
    }
    if (any(ypos2)) {
      h@entries[ypos2] <- theta[2] * (1 - theta[3]) * ((1 - tmp[ypos2])^6 * 
                                        (35 * tmp[ypos2]^2 + 18 * tmp[ypos2] + 3))/3
    }
    if (any(npos)) {
      h@entries[npos] <- 0
    }
  }
  else {
    tmp <- c(h)/theta[1]
    ypos <- tmp < eps
    ypos2 <- tmp >= eps & tmp < 1
    npos <- tmp >= 1
    if (any(ypos)) {
      h[ypos] <- theta[2]
    }
    if (any(ypos2)) {
      h[ypos2] <- theta[2] * (1 - theta[3]) * ((1 - tmp[ypos2])^6 * (35 * 
                                                      tmp[ypos2]^2 + 18 * tmp[ypos2] + 3))/3
    }
    if (any(npos)) {
      h[npos] <- 0
    }
  }
  return(h)
}