# Figure 1 ----

{
  rm(list = ls())
  
  source('00_system.R')
  
  colors_to_use <- palette.colors(palette = 'Okabe-Ito')
  do_analysis <- F
  
  C_de_y <- function(z_transf, kappa, alpha){
    cosh(kappa * asinh(z_transf) - alpha)
  }
  
  S_de_y <- function(z_transf, alpha, kappa){
    sinh(kappa * asinh(z_transf) - alpha)
  }
  
  density_SAS <- function(x, kappa, alpha){
    kappa * C_de_y(x, kappa = kappa, alpha = alpha) / sqrt(2 * pi * (1 + x^2)) * 
      exp(-S_de_y(z_transf = x, alpha = alpha, kappa = kappa)^2/2)
  }
  
  curve(density_SAS(x, kappa = 1.2, alpha = 0), from = -5, to = 5, n = 1000)
  
  side_left <- c(0.4, 0.8, 1, 1.2)
  
  x_seq <- seq(-10, 10, length.out = 1000)
  
  colors_to_use <- palette.colors(palette = 'Okabe-Ito')
  
  {
    tikz(paste0(getwd(), '/finalplots/SAS_density_left.tex'), standAlone = T, width = 7.5, height = 3.75)
    par(pty = 'm', oma = c(0, 3, 0, 0), mar = c(2, 0, 1.5, 0.5), mfrow = c(1, 1))
    
    for(i in 1:length(side_left)){
      
      if(i == 1){
        plot(x = x_seq, y = density_SAS(x_seq, kappa = side_left[i], alpha = 0),
             type = 'l' , ylab = ' ', col = colors_to_use[i], ylim = c(0, 0.5), xlab = '', lwd = 7, axes = F)
      } else {
        lines(x = x_seq, y = density_SAS(x_seq, kappa = side_left[i], alpha = 0),
              ylab = ' ', col = colors_to_use[i], lwd = 7)
      }
    }
    
    box()
    axis(1, cex.axis = 2)
    axis(2, cex.axis = 2)
    
    legend(bty = 'n', 'topright', legend = c(paste0('$\\alpha = 0$', ', $\\kappa$ = ', side_left[1]),
                                             paste0('$\\alpha = 0$', ', $\\kappa$ = ', side_left[2]),
                                             paste0('$\\alpha = 0$', ', $\\kappa$ =  ', side_left[3]),
                                             paste0('$\\alpha = 0$', ', $\\kappa$ = ', side_left[4]))
           , col = (colors_to_use)[1:4], lty = c(1, 1, 1, 1), cex = 1.5)
    
    dev.off()
    tinytex::pdflatex(paste0(CommonPath, '/finalplots/SAS_density_left.tex'))
    system(paste0('rm ~', '/finalplots/SAS_density_left.tex'))
  }
  
  {
    right_side <- c(0.4, 0.8, 1, 1.2)
    
    tikz(paste0(CommonPath, '/finalplots/SAS_density_right.tex'), standAlone = T, width = 7.5, height = 3.75)
    par(pty = 'm', oma = c(0, 3, 0, 0), mar = c(2, 0, 1.5, 0.5), mfrow = c(1, 1))
    
    
    for(i in 1:length(right_side)){
      
      if(i == 1){
        plot(x = x_seq, y = density_SAS(x_seq, kappa = 1, alpha = right_side[i]),
             type = 'l' , ylab = ' ', col = colors_to_use[i], ylim = c(0, 0.4), xlim = c(-4, 8), xlab = '', lwd = 7, axes = F)
      } else {
        lines(x = x_seq, y = density_SAS(x_seq, kappa = 1, alpha = right_side[i]),
              ylab = ' ', col = colors_to_use[i], lwd = 7)
      }
    }
    
    box()
    axis(1, cex.axis = 2)
    axis(2, cex.axis = 2)
    legend(bty = 'n', 'topright', legend = c(paste0('$\\alpha =$ ', right_side[1], ', $\\kappa$ = 1'),
                                             paste0('$\\alpha =$ ', right_side[2], ', $\\kappa$ = 1'),
                                             paste0('$\\alpha =$ ', '1.0', ', $\\kappa$ = 1'),
                                             paste0('$\\alpha =$ ', right_side[4], ', $\\kappa$ = 1'))
           , col = (colors_to_use)[1:4], lty = c(1, 1, 1, 1), cex = 1.5)
    
    dev.off()
    
    tinytex::pdflatex(paste0(CommonPath, '/finalplots/SAS_density_right.tex'))
    system(paste0('rm ~', '/finalplots/SAS_density_right.tex'))
  }
  
  rm(list = ls())
  
}

# Figures from simulation study setup ----

rm(list = ls())
source('00_system.R')

colors_to_use <- palette.colors(palette = 'Okabe-Ito')
do_analysis <- F

list_estimates <- list.files(paste0(CommonPath, '/estimates'))
list_estimates <- list_estimates[grep('estimates.', list_estimates)]

to_fill <- matrix(NA, ncol = 4, nrow = 200)
colnames(to_fill) <- c('Matern', 'Tapering', 'Directmiss', 'Composite')

# Thetas ----

result_theta_1 <- lapply(1:16, function(x) to_fill)
result_theta_2 <- lapply(1:16, function(x) to_fill)
result_theta_4 <- lapply(1:16, function(x) to_fill)

for(i in 1:length(list_estimates)){
  
  load(paste0(CommonPath, '/estimates/', list_estimates[i]))
  load(paste0(CommonPath, '/fields/', this_estimates@loadstring))
  
  if(is(this_estimates, 'directmiss_estimate')){
    
    z <- 3
    
    result_theta_2[[this_estimates@specsnr]][this_estimates@fielditer, z] <-  this_estimates@optimout$par[1]
    result_theta_4[[this_estimates@specsnr]][this_estimates@fielditer, z] <-  this_estimates@optimout$par[2]
    
  } else{
    
    if(is(this_estimates, 'composite_estimate')){
      
      z <- 4
      
      if(this_field@gaussian){
        result_theta_1[[this_estimates@specsnr]][this_estimates@fielditer, z] <-  unlist(this_estimates@optimout$par[2])
        result_theta_2[[this_estimates@specsnr]][this_estimates@fielditer, z] <-  unlist(this_estimates@optimout$par[3])
        result_theta_4[[this_estimates@specsnr]][this_estimates@fielditer, z] <-  unlist(this_estimates@optimout$par[1])
        
      } else{
        result_theta_1[[this_estimates@specsnr]][this_estimates@fielditer, z] <-  this_estimates@optimout$par[5]
        result_theta_2[[this_estimates@specsnr]][this_estimates@fielditer, z] <-  this_estimates@optimout$par[6]
        result_theta_4[[this_estimates@specsnr]][this_estimates@fielditer, z] <-  this_estimates@optimout$par[4]
      }
      
    } else{
      
      if(is(this_estimates, 'taper_estimate')){z <- 2}
      if(is(this_estimates, 'matern_estimate')){z <- 1}
      
      result_theta_1[[this_estimates@specsnr]][this_estimates@fielditer, z] <-  this_estimates@optimout$par[1]
      result_theta_2[[this_estimates@specsnr]][this_estimates@fielditer, z] <-  this_estimates@optimout$par[2]
      result_theta_4[[this_estimates@specsnr]][this_estimates@fielditer, z] <-  this_estimates@optimout$par[3]
      
    }
    
  }
  
  if(i == 1){
    
    bias_analysis_theta_1 <- data.frame('bias' = (as.numeric(result_theta_1[[this_estimates@specsnr]][this_estimates@fielditer, z]-this_field@truetheta[1])/this_field@truetheta[1]),
                                        'gridded' = this_field@gridded, 'gauss' = this_field@gaussian, 'size' = factor(as.character(this_field@nsize), levels = c('784', '3025')),
                                        'thetas' = factor(this_field@parsetting, levels = c('low', 'high')), 'estim' = factor(as.character(z), levels = c('1', '2', '3', '4')))
    
    bias_analysis_theta_2 <- data.frame('bias' = (as.numeric(result_theta_2[[this_estimates@specsnr]][this_estimates@fielditer, z]-this_field@truetheta[2])/this_field@truetheta[2]),
                                        'gridded' = this_field@gridded, 'gauss' = this_field@gaussian, 'size' = factor(as.character(this_field@nsize), levels = c('784', '3025')),
                                        'thetas' = factor(this_field@parsetting, levels = c('low', 'high')), 'estim' = factor(as.character(z), levels = c('1', '2', '3', '4')))
    
    bias_analysis_theta_4 <- data.frame('bias' = (as.numeric(result_theta_4[[this_estimates@specsnr]][this_estimates@fielditer, z]-(this_field@truetheta[4]))),
                                        'gridded' = this_field@gridded, 'gauss' = this_field@gaussian, 'size' = factor(as.character(this_field@nsize), levels = c('784', '3025')),
                                        'thetas' = factor(this_field@parsetting, levels = c('low', 'high')), 'estim' = factor(as.character(z), levels = c('1', '2', '3', '4')))
    
  } else {
    
    bias_analysis_theta_1 <- rbind.data.frame(bias_analysis_theta_1, setNames(
      data.frame((as.numeric(result_theta_1[[this_estimates@specsnr]][this_estimates@fielditer, z]-this_field@truetheta[1])/this_field@truetheta[1]),
                 this_field@gridded,
                 this_field@gaussian,
                 this_field@nsize,
                 this_field@parsetting,
                 z), names(bias_analysis_theta_1)))
    
    bias_analysis_theta_2 <- rbind.data.frame(bias_analysis_theta_2, setNames(
      data.frame((as.numeric(result_theta_2[[this_estimates@specsnr]][this_estimates@fielditer, z]-this_field@truetheta[2])/this_field@truetheta[2]),
                 this_field@gridded,
                 this_field@gaussian,
                 this_field@nsize,
                 this_field@parsetting,
                 z), names(bias_analysis_theta_2)))
    
    bias_analysis_theta_4 <- rbind.data.frame(bias_analysis_theta_4, setNames(
      data.frame((as.numeric((result_theta_4[[this_estimates@specsnr]][this_estimates@fielditer, z]-(this_field@truetheta[4])))),
                 this_field@gridded,
                 this_field@gaussian,
                 this_field@nsize,
                 this_field@parsetting,
                 z), names(bias_analysis_theta_4)))
    
  }
  
}

# Betas ----

result_beta_0 <- lapply(1:16, function(x) to_fill)
bias_analysis_beta_0 <- data.frame('bias' = NA_real_, 'gridded' = NA, 'gauss' = NA, 'size' = NA_integer_, 'thetas' = NA, 'estim' = NA)

for(i in 1:length(list_estimates)){
  load(paste0(CommonPath, '/estimates/', list_estimates[i]))
  load(paste0(CommonPath, '/fields/', this_estimates@loadstring))
  
  if(is(this_estimates, 'composite_estimate')){z <- 4}
  if(is(this_estimates, 'directmiss_estimate')){z <- 3}
  if(is(this_estimates, 'taper_estimate')){z <- 2}
  if(is(this_estimates, 'matern_estimate')){z <- 1}
  
  result_beta_0[[this_estimates@specsnr]][this_estimates@fielditer, z] <- this_estimates@betas_hat[1]

  if(i == 1){
    
    bias_analysis_beta_0 <- data.frame('bias' = (as.numeric(result_beta_0[[this_estimates@specsnr]][this_estimates@fielditer, z]-this_field@truebetas[1])/this_field@truebetas[1]),
                                       'gridded' = this_field@gridded, 'gauss' = this_field@gaussian, 'size' = factor(as.character(this_field@nsize), levels = c('784', '3025')),
                                       'thetas' = factor(this_field@parsetting, levels = c('low', 'high')), 'estim' = factor(as.character(z), levels = c('1', '2', '3', '4')))
    
  } else {
    
    bias_analysis_beta_0 <- rbind.data.frame(bias_analysis_beta_0, setNames(
      data.frame((as.numeric(result_beta_0[[this_estimates@specsnr]][this_estimates@fielditer, z]-this_field@truebetas[1])/this_field@truebetas[1]),
                 this_field@gridded,
                 this_field@gaussian,
                 this_field@nsize,
                 this_field@parsetting,
                 z), names(bias_analysis_beta_0)))
    
  }
}

# Computational measures ----

result_memory <- lapply(1:16, function(x) to_fill)
result_time <- lapply(1:16, function(x) to_fill)

for(i in 1:length(list_estimates)){
  load(paste0(CommonPath, '/estimates/', list_estimates[i]))
  load(paste0(CommonPath, '/fields/', this_estimates@loadstring))
  
  if(is(this_estimates, 'composite_estimate')){z <- 4}
  if(is(this_estimates, 'directmiss_estimate')){z <- 3}
  if(is(this_estimates, 'taper_estimate')){z <- 2}
  if(is(this_estimates, 'matern_estimate')){z <- 1}
  
  result_memory[[this_estimates@specsnr]][this_estimates@fielditer, z] <-  this_estimates@memory_used
  result_time[[this_estimates@specsnr]][this_estimates@fielditer, z] <-  this_estimates@time
  
}


# Figure 2 - betas 0----

if(do_analysis){
  
  plot(density(bias_analysis_beta_0$bias[complete.cases(bias_analysis_beta_0$bias)]))
  model_beta_0 <- lm(bias ~ . + estim*thetas + estim*gridded + estim*gauss + estim*size, data = bias_analysis_beta_0)
  anova(model_beta_0)
  plot(model_beta_0)
}

{
  
  mean_values_beta_0 <- aggregate(bias~estim, data = bias_analysis_beta_0, FUN = mean)$'bias'
  mean_values_beta_0_b <- aggregate(bias~estim + thetas, data = bias_analysis_beta_0, FUN = mean)$'bias'
  boxplot(data = bias_analysis_beta_0, bias ~  gauss + size + gridded +  thetas + estim)
  
  text_names <- c('beta0_a.tex', 'beta0_b.tex', 'beta0_c.tex', 'beta0_d.tex')
  
  lower_bound <- -5.5
  
  for(i in 1:4){
    
    tikz(paste0(CommonPath, '/finalplots/', text_names[i]), standAlone = T, width = 4, height = 6)
    
    par(pty = 'm', oma = c(3, 3, 0, 0), mar = c(2, 0, 1.5, 0.5), mfrow = c(1, 1))
    plot(0, type = 'n', axes = FALSE, ann = FALSE, ylim= c(lower_bound, 6), xlim = c(0.5, 4.5))
    axis(side = 2, cex.axis = 2, tck = 1, lty = 2, col = 'grey70', labels = ifelse(i == 1, T, F))
    box()
    
    boxplot(data = bias_analysis_beta_0[which(bias_analysis_beta_0$estim == i), ], bias ~  gauss + thetas, 
            axes = FALSE, cex.box = 2, boxlwd = 4, outwex = 3, boxwex = 0.75, add = T, main = '')
    
    if(i == 1) {
      segments(0.5, mean_values_beta_0[1], 4.5, mean_values_beta_0[1], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_beta_0_b[1], 2.5, mean_values_beta_0_b[1], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_beta_0_b[5], 4.5, mean_values_beta_0_b[5], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    if(i == 2){
      segments(0.5, mean_values_beta_0[2], 4.5, mean_values_beta_0[2], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_beta_0_b[2], 2.5, mean_values_beta_0_b[2], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_beta_0_b[6], 4.5, mean_values_beta_0_b[6], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    if(i == 3){
      segments(0.5, mean_values_beta_0[2], 4.5, mean_values_beta_0[2], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_beta_0_b[3], 2.5, mean_values_beta_0_b[3], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_beta_0_b[7], 4.5, mean_values_beta_0_b[7], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    if(i == 4){
      segments(0.5, mean_values_beta_0[3], 4.5, mean_values_beta_0[3], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_beta_0_b[4], 2.5, mean_values_beta_0_b[4], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_beta_0_b[8], 4.5, mean_values_beta_0_b[8], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    
    # first line labels
    axis(side = 1, at = c(1, 2, 3, 4), cex.axis = 2, labels = c('SAS', 'Gau', 'SAS', 'Gau'))
    
    # second line labels
    
    axis(side = 1, at = c(1, 1.5, 2), cex.axis = 2, labels = c(NA, NA, NA), outer = T, tck = 0, pos = lower_bound - 1.6)
    
    axis(side = 1, at = c(1.5), cex.axis = 2, labels = c('Low'), outer = T, tck = F, pos = lower_bound - 1.6)
    
    axis(side = 1, at = c(1.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos = lower_bound - 1.6)
    axis(side = 1, at = c(1, 2), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound - 1.6)
    
    axis(side = 1, at = c(3, 3.5, 4), cex.axis = 2, labels = c(NA, NA, NA), outer = T, tck = 0, pos = lower_bound - 1.6)
    axis(side = 1, at = c(3.5), cex.axis = 2, labels = c('High'), outer = T, tck = F, pos = lower_bound - 1.6)
    axis(side = 1, at = c(3.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos = lower_bound - 1.6)
    axis(side = 1, at = c(3,4), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound - 1.6)
    
    
    dev.off()
    
    tinytex::pdflatex(paste0(CommonPath, '/finalplots/', text_names[i]))
    
    system(paste0('rm ~', '/finalplots/', text_names[i]))
  
  }
}

# Figure 3 - scale parameter ----

if(do_analysis){
  
  plot(density(bias_analysis_theta_2$bias))
  model_theta_2 <- lm(bias ~ . + estim*thetas + estim*gridded + estim*gauss + estim*size, data = bias_analysis_theta_2)
  anova(model_theta_2)
  model_theta_2
  plot(model_theta_2)
  
}

{
  mean_values_theta_2 <- aggregate(bias~estim, data = bias_analysis_theta_2, FUN = mean)$'bias'
  mean_values_theta_2_b <- aggregate(bias~estim*thetas, data = bias_analysis_theta_2, FUN = mean)$'bias'
  boxplot(data = bias_analysis_theta_2, bias ~ gauss + size + gridded + thetas + estim, main = '  ')
  
  text_names <- c('theta2_a.tex', 'theta2_b.tex', 'theta2_c.tex', 'theta2_d.tex')
  
  lower_bound <- -1
  
  for(i in 1:4){
    
    tikz(paste0(CommonPath, '/finalplots/', text_names[i]), standAlone = T, width = 4, height = 6)
    
    par(pty = 'm', oma = c(3, 3, 0, 0), mar = c(2, 0, 1.5, 0.5), mfrow = c(1, 1))
    
    plot(0, type = 'n', axes = FALSE, ann = FALSE, ylim = c(lower_bound, 7.5), xlim = c(0.5, 4.5))
    axis(side = 2, cex.axis = 2, tck = 1, lty = 2, col = 'grey70', labels = ifelse(i == 1, T, F))
    box()
    
    
    boxplot(data = bias_analysis_theta_2[which(bias_analysis_theta_2$estim == i), ], bias ~ gauss + thetas, main = ' ',
            axes = FALSE, cex.box = 2, boxlwd = 4, outwex = 3, boxwex = 0.75, add = T)
    
    if(i == 1){
      segments(0.5, mean_values_theta_2[1], 4.5, mean_values_theta_2[1], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_theta_2_b[1], 2.5, mean_values_theta_2_b[1], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_theta_2_b[5], 4.5, mean_values_theta_2_b[5], col = colors_to_use[6], lty = 2, lwd = 4)
      
    }
    
    if(i == 2){
      segments(0.5, mean_values_theta_2[2], 4.5, mean_values_theta_2[2], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_theta_2_b[2], 2.5, mean_values_theta_2_b[2], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_theta_2_b[6], 4.5, mean_values_theta_2_b[6], col = colors_to_use[6], lty = 2, lwd = 4)
      
    }
    
    if(i == 3){
      segments(0.5, mean_values_theta_2[3], 4.5, mean_values_theta_2[3], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_theta_2_b[3], 2.5, mean_values_theta_2_b[3], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_theta_2_b[7], 4.5, mean_values_theta_2_b[7], col = colors_to_use[6], lty = 2, lwd = 4)
      
    }
    
    if(i == 4){
      segments(0.5, mean_values_theta_2[4], 4.5, mean_values_theta_2[4], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_theta_2_b[4], 2.5, mean_values_theta_2_b[4], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_theta_2_b[8], 4.5, mean_values_theta_2_b[8], col = colors_to_use[6], lty = 2, lwd = 4)
      
    }
    
    # first line labels
    axis(side = 1, at = c(1, 2, 3, 4), cex.axis = 2, labels = c('SAS', 'Gau', 'SAS', 'Gau'))
    
    # second line labels
    
    axis(side = 1, at = c(1, 1.5, 2), cex.axis = 2, labels = c(NA, NA, NA), outer = T,tck = 0, pos = lower_bound - 1.1)
    
    axis(side = 1, at = c(1.5), cex.axis = 2, labels = c('Low'), outer = T, tck = F, pos = lower_bound - 1.1)
    
    axis(side = 1, at = c(1.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos = lower_bound - 1.1)
    axis(side = 1, at = c(1,2), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound - 1.1)
    
    axis(side = 1, at = c(3, 3.5, 4), cex.axis = 2, labels = c(NA, NA, NA), outer = T, tck = 0, pos = lower_bound - 1.1)
    axis(side = 1, at = c(3.5), cex.axis = 2, labels = c('High'), outer = T, tck = F, pos = lower_bound - 1.1)
    axis(side = 1, at = c(3.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos = lower_bound - 1.1)
    axis(side = 1, at = c(3, 4), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound - 1.1)
    
    dev.off()
    
    tinytex::pdflatex(paste0(CommonPath, '/finalplots/', text_names[i]))
    system(paste0('rm ~', '/finalplots/', text_names[i]))
    
    
  }
  
}

# Figure 4 - range -----

if(do_analysis){
  
  plot(density(bias_analysis_theta_1$bias[complete.cases(bias_analysis_theta_1$bias)]))
  plot(density(sqrt(abs(bias_analysis_theta_1$bias[complete.cases(bias_analysis_theta_1$bias)]/this_field@truetheta[1]))))
  model_theta_1 <- lm(sqrt(abs(bias/this_field@truetheta[1])) ~ . + estim*thetas + estim*gridded + estim*gauss + estim*size, data = bias_analysis_theta_1)
  anova(model_theta_1)
  model_theta_1
  plot(model_theta_1)
  
}

{
  
  mean_values_theta_1 <- aggregate(sqrt(abs(bias/this_field@truetheta[1]))~estim + thetas, data = bias_analysis_theta_1, FUN = mean)$'sqrt(abs(bias/this_field@truetheta[1]))'
  mean_values_theta_1_b <- aggregate(sqrt(abs(bias/this_field@truetheta[1]))~ estim, data = bias_analysis_theta_1, FUN = mean)$'sqrt(abs(bias/this_field@truetheta[1]))'
  boxplot(data = bias_analysis_theta_1, sqrt(abs(bias/this_field@truetheta[1])) ~ gauss + size + gridded + thetas + estim)
  
  # Range figure
  text_names <- c('theta1_a.tex', 'theta1_b.tex', '', 'theta1_c.tex')
  num_index <- c( 1, 2, 4)
  lower_bound <- -0.1
  for(i in num_index){
    
    tikz(paste0(CommonPath, '/finalplots/', text_names[i]), standAlone = T, width = 4, height = 6)
    
    par(pty = 'm', oma = c(3, 3, 0, 0), mar = c(2, 0, 1.5, 0.5), mfrow = c(1, 1))
    if(i == 1 | i == 4){
      plot(0, type = 'n', axes = FALSE, ann = FALSE, ylim = c(lower_bound, 30), xlim = c(0.5, 4.5))
      
    } else (plot(0, type = 'n', axes = FALSE, ann = FALSE, ylim = c(lower_bound, 2600), xlim = c(0.5, 4.5)))
    axis(side = 2, cex.axis = 2, tck = 1, lty = 2, col = 'grey70', labels = T)
    box()
    
    
    boxplot(data = bias_analysis_theta_1[which(bias_analysis_theta_1$estim == i), ], sqrt(abs(bias/this_field@truetheta[1])) ~ size + thetas, main = ' ',
            axes = FALSE, cex.box = 2, boxlwd = 4, outwex = 3, boxwex = 0.75, add = T)
    
    if(i == 1 | i == 4){
      # first line labels
      axis(side = 1, at = c(1, 2, 3, 4), cex.axis = 2, labels = c('784', '3025', '784', '3025'))
      
      # second line labels
      
      axis(side = 1, at = c(1, 1.5, 2), cex.axis = 2, labels = c(NA, NA, NA), outer = T, tck = 0, pos = lower_bound - 4)
      
      axis(side = 1, at = c(1.5), cex.axis = 2, labels = c('Low'), outer = T, tck = F, pos = lower_bound - 4)
      
      axis(side = 1, at = c(1.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos = lower_bound - 4)
      axis(side = 1, at = c(1, 2), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound - 4)
      
      axis(side = 1, at = c(3, 3.5, 4), cex.axis = 2, labels = c(NA, NA, NA), outer = T, tck = 0, pos = lower_bound - 4)
      
      axis(side = 1, at = c(3.5), cex.axis = 2, labels = c('High'), outer = T, tck = F, pos = lower_bound - 4)
      axis(side = 1, at = c(3.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos =lower_bound - 4)
      axis(side = 1, at = c(3, 4), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound - 4)
      
    } else {
      
      axis(side = 1, at = c(1, 2, 3, 4), cex.axis = 2, labels = c('784', '3025', '784', '3025'))
      
      # second line labels
      
      axis(side = 1, at = c(1, 1.5, 2), cex.axis = 2, labels = c(NA, NA, NA), outer = T, tck = 0, pos = lower_bound -5.5*60)
      
      axis(side = 1, at = c(1.5), cex.axis = 2, labels = c('Low'), outer = T, tck = F, pos = lower_bound - 5.5*60)
      
      axis(side = 1, at = c(1.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos = lower_bound -5.5*60)
      axis(side = 1, at = c(1, 2), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound -5.5*60)
      
      axis(side = 1, at = c(3, 3.5, 4), cex.axis = 2, labels = c(NA, NA, NA), outer = T, tck = 0, pos = lower_bound -5.5*60)
      
      axis(side = 1, at = c(3.5), cex.axis = 2, labels = c('High'), outer = T, tck = F, pos = lower_bound -5.5*60)
      axis(side = 1, at = c(3.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos = lower_bound -5.5*60)
      axis(side = 1, at = c(3, 4), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound -5.5*60)
      
    }
    
    if(i == 1) {
      segments(0.5, mean_values_theta_1_b[1], 4.5, mean_values_theta_1_b[1], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_theta_1[1], 2.5, mean_values_theta_1[1], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_theta_1[4], 4.5, mean_values_theta_1[4], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    if(i == 2){
      segments(0.5, mean_values_theta_1_b[2], 4.5, mean_values_theta_1_b[2], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_theta_1[2], 2.5, mean_values_theta_1[2], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_theta_1[5], 4.5, mean_values_theta_1[5], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    if(i == 4){
      segments(0.5, mean_values_theta_1_b[3], 4.5, mean_values_theta_1_b[3], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_theta_1[3], 2.5, mean_values_theta_1[3], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_theta_1[6], 4.5, mean_values_theta_1[6], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    
    
    dev.off()
    
    tinytex::pdflatex(paste0(CommonPath, '/finalplots/', text_names[i]))
    system(paste0('rm ~', '/finalplots/', text_names[i]))
  }
  
}

# Figure 5 - white-noise -----

if(do_analysis){
  
  plot(density(bias_analysis_theta_4$bias))
  model_theta_4 <- lm(bias ~ . + estim*thetas + estim*gridded + estim*gauss, data = bias_analysis_theta_4)
  anova(model_theta_4)
  model_theta_4
  plot(model_theta_4)
  
}

{
  
  mean_values_theta_4 <- aggregate(bias~estim, data = bias_analysis_theta_4, FUN = mean)$'bias'
  mean_values_theta_4_b <- aggregate(bias~ estim + thetas, data = bias_analysis_theta_4, FUN = mean)$'bias'
  
  text_names <- c('theta4_a.tex', 'theta4_b.tex', 'theta4_c.tex', 'theta4_d.tex')
  lower_bound <- -0.1
  for(i in 1:4){
    
    tikz(paste0(CommonPath, '/finalplots/', text_names[i]), standAlone = T, width = 4, height = 6)
    
    par(pty = 'm', oma = c(3, 3, 0, 0), mar = c(2, 0, 1.5, 0.5), mfrow = c(1, 1))
    
    plot(0, type = 'n', axes = FALSE, ann = FALSE, ylim = c(lower_bound, 0.6), xlim = c(0.5, 4.5))
    axis(side = 2, cex.axis = 2, tck = 1, lty = 2, col = 'grey70', labels = ifelse(i == 1, T, F))
    box()
    
    
    boxplot(data = bias_analysis_theta_4[which(bias_analysis_theta_4$estim == i), ], bias ~ size + thetas, main = ' ',
            axes = FALSE, cex.box = 2, boxlwd = 4, outwex = 3, boxwex = 0.75, add = T)
    
    if(i == 1){
      segments(0.5, mean_values_theta_4[1], 4.5, mean_values_theta_4[1], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_theta_4_b[1], 2.5, mean_values_theta_4_b[1], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_theta_4_b[5], 4.5, mean_values_theta_4_b[5], col = colors_to_use[6], lty = 2, lwd = 4)
      
    }
    
    if(i == 2){
      segments(0.5, mean_values_theta_4[2], 4.5, mean_values_theta_4[2], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_theta_4_b[2], 2.5, mean_values_theta_4_b[2], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_theta_4_b[6], 4.5, mean_values_theta_4_b[6], col = colors_to_use[6], lty = 2, lwd = 4)
      
    }
    
    if(i == 3){
      segments(0.5, mean_values_theta_4[3], 4.5, mean_values_theta_4[3], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_theta_4_b[3], 2.5, mean_values_theta_4_b[3], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_theta_4_b[7], 4.5, mean_values_theta_4_b[7], col = colors_to_use[6], lty = 2, lwd = 4)
      
    }
    
    if(i == 4){
      segments(0.5, mean_values_theta_4[4], 4.5, mean_values_theta_4[4], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_theta_4_b[4], 2.5, mean_values_theta_4_b[4], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_theta_4_b[8], 4.5, mean_values_theta_4_b[8], col = colors_to_use[6], lty = 2, lwd = 4)
      
    }
    
    # first line labels
    axis(side = 1, at = c(1, 2, 3, 4), cex.axis = 2, labels = c('784', '3025', '784', '3025'))
    
    # second line labels
    
    axis(side = 1, at = c(1, 1.5, 2), cex.axis = 2, labels = c(NA, NA, NA), outer = T, tck = 0, pos = lower_bound -0.09)
    
    axis(side = 1, at = c(1.5), cex.axis = 2, labels = c('Low'), outer = T, tck = F, pos = lower_bound -0.09)
    
    axis(side = 1, at = c(1.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos = lower_bound -0.09)
    axis(side = 1, at = c(1, 2), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound -0.09)
    
    axis(side = 1, at = c(3, 3.5, 4), cex.axis = 2, labels = c(NA, NA, NA), outer = T, tck = 0, pos = lower_bound -0.09)
    axis(side = 1, at = c(3.5), cex.axis = 2, labels = c('High'), outer = T, tck = F, pos = lower_bound -0.09)
    axis(side = 1, at = c(3.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos = lower_bound -0.09)
    axis(side = 1, at = c(3, 4), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound -0.09)
    
    dev.off()
    
    tinytex::pdflatex(paste0(CommonPath, '/finalplots/', text_names[i]))
    system(paste0('rm ~', '/finalplots/', text_names[i]))

  }
}

# Figure 6 - computational summary ----

tmp_time <- result_time[[1]]
tmp_time <- tmp_time/tmp_time[, 1]
time_boxplot <- cbind(tmp_time, rep(1, dim(tmp_time)[1]))
for(i in 2:16){ 
  tmp_time <- result_time[[i]]
  tmp_time <- tmp_time/tmp_time[, 1]
  if(i %in% c(1, 2, 3, 4, 9, 10, 11, 12)) {
    tmp_time <- cbind(tmp_time, rep(1, dim(tmp_time)[1]))}
  else tmp_time <- cbind(tmp_time, rep(2, dim(tmp_time)[1]))
  
  time_boxplot <- rbind(time_boxplot, tmp_time)
}

sum(time_boxplot[, 2] > 1) / dim(time_boxplot)[1] * 100
sum(time_boxplot[, 3] > 1) / dim(time_boxplot)[1] * 100
sum(time_boxplot[, 4] > 1) / dim(time_boxplot)[1] * 100

{
  tikz(paste0(CommonPath, '/finalplots/', 'time_percentage.tex'), standAlone = T, width = 4, height = 6)
  par(pty = 'm', oma = c(0, 0, 0, 0), mar = c(2, 2.5, 1.5, 2.5))
  plot(0:1, 0:1, type = 'n', xlim = c( 0.5, 3.5), ylim = c( 0, 3),
       axes = FALSE, ann = FALSE)
  axis(side = 2, cex.axis = 2, tck = 1, lty = 2, col = 'grey70')
  box()
  axis(side = 1, at = c(1, 2, 3), cex.axis = 2, labels = c('TA', 'DM', 'CL'))
  
  vioplot((time_boxplot[which(time_boxplot[, 5] == 1 ), c(2:4)]), h = 0.03, areaEqual = F, wex = 1.2, col = 'gray70',
          add = T, side = 'left', axes = F, ylab = NA)
  vioplot((time_boxplot[which(time_boxplot[, 5] == 2), 2:4]), h = 0.03, areaEqual = F, wex = 1.2, col = 'gray70',
          add = T, side = 'right', axes = F)
  abline(h = 1, col = colors_to_use[3], lty = 2 , lwd = 5)
  dev.off()
  tinytex::pdflatex(paste0(CommonPath, '/finalplots/', 'time_percentage.tex'))
  system(paste0('rm ~', '/finalplots/', 'time_percentage.tex'))

}

tmp_memory <- result_memory[[1]]
tmp_memory <- tmp_memory/tmp_memory[, 1]
memory_boxplot <- cbind(tmp_memory, rep(1, dim(tmp_memory)[1]))
for(i in 2:16){
  tmp_memory <- result_memory[[i]]
  tmp_memory <- tmp_memory/tmp_memory[, 1]
  if(i %in% c(1, 2, 3, 4, 9, 10, 11, 12)) {
    tmp_memory <- cbind(tmp_memory, rep(1, dim(tmp_memory)[1]))}   else tmp_memory <- cbind(tmp_memory, rep(2, dim(tmp_memory)[1]))
  
  memory_boxplot <- rbind(memory_boxplot, tmp_memory)
}

{
  tikz(paste0(CommonPath, '/finalplots/', 'memory_percentage.tex'), standAlone = T, width = 4, height = 6)
  
  par(pty = 'm', oma = c(0, 0, 0, 0), mar = c(2, 2.5, 1.5, 2.5))
  plot(0:1, 0:1, type = 'n', xlim = c(0.5, 3.5), ylim = c(0.1, 2),
       axes = FALSE, ann = FALSE)
  axis(side = 2, at = c(0.2, 0.6, 1, 1.6, 2), cex.axis = 2, tck = 1, lty = 2, col = 'grey70')
  box()
  axis(side = 1, at = c(1, 2, 3), cex.axis = 2, labels = c('TA', 'DM', 'CL'))
  vioplot(memory_boxplot[which(memory_boxplot[, 5] == 1), 2:4], h = 0.04, areaEqual = F, wex = 1.25, ylim = c(-6, 3), col = 'gray70', add = T, side = 'left', axes = F, ylab = NA)
  vioplot(memory_boxplot[which(memory_boxplot[, 5] == 2), 2:4], h=0.04, areaEqual = F, wex = 1.25, ylim = c(-6, 3), col = 'gray70', add = T, side = 'right', axes  =  F)
  abline(h = 1, col = colors_to_use[3], lty = 2 , lwd = 5)
  dev.off()
  
  tinytex::pdflatex(paste0(CommonPath, '/finalplots/', 'memory_percentage.tex'))
  system(paste0('rm ~', '/finalplots/', 'memory_percentage.tex'))

}

# Git figures

# Figure 7 - maps ----

rm(list = ls())
source('00_system.R')
load(paste0(CommonPath, '/application/ClimChange.RData'))

{
  pdf(paste0(CommonPath, '/finalplots/illustration_left.pdf'), width = 8, height = 4*1.25)
  par(pty = 'm', oma = c(0, 0, 0, 0), mar = c(1, 0, 0, 0), mfrow = c(1, 1))
  image.plot(long, lati, tas1, axes = F, ylab = '', xlab = '', axis.args=list(cex.axis = 1.5), horizontal = T, legend.width = 0.5)
  box() 
  map(add = TRUE, wrap = c(0, 360))
  dev.off()
}

{
  pdf(paste0(CommonPath, '/finalplots/illustration_right.pdf'), width = 8, height = 4*1.25)
  par(pty = 'm', oma = c(0, 0, 0, 0), mar = c(1, 0, 0, 0), mfrow = c(1, 1))
  image.plot(long, lati, tas2 - tas1, axes = F, ylab = '', xlab = '', axis.args = list(cex.axis = 1.5), horizontal = T, legend.width = 0.5)
  box()
  map(add = TRUE, wrap = c(0, 360))
  dev.off()
}

# Figure 8 - illustration results ----

rm(list = ls())
source('00_system.R')
colors_to_use <- palette.colors(palette = 'Okabe-Ito')

load(file = paste0(CommonPath, '/application/illustration.RData'))

{
  z_values <- rnorm(1000000, mean = 0, sd = 1)
  tikz(paste0(CommonPath, '/finalplots/illustration_results.tex'), standAlone = T, width = 7.5, height = 3.75)
  par(pty = 'm', oma = c(0, 2.5, 0, 0), mar = c(2, 0, 1.5, 1), mfrow = c(1, 1))
  plot(density(c(y)), axes = FALSE, ann = FALSE, xlim = c(-5, 8), ylim = c(0, 0.5), lwd = 7)
  box()
  axis(1, cex.axis = 2)
  axis(2, cex.axis = 2)
  lines(density(sqrt(output_geofit$param[3])*(sinh((1/output_geofit$param[5])*(asinh(z_values) - output_geofit$param[4])))), 
        col = colors_to_use[3], lty = 1, lwd = 7)
  lines(density(sqrt(output_taper$par[2])*(sinh((1/output_taper$par[5])*(asinh(z_values) - output_taper$par[4])))),
        col = colors_to_use[2], lty = 3, lwd = 7)
  legend('topright', bty = 'n', legend = c('Data', 'CL', 'TA $\\&$ DM'), col = c(colors_to_use[1], colors_to_use[6], col = colors_to_use[2]), 
         lty = c(1, 1, 3), lwd = c(7, 7, 7), cex = 1.75)
  dev.off()
  tinytex::pdflatex(paste0(CommonPath, '/finalplots/illustration_results.tex'))
  system(paste0('rm ~', '/finalplots/illustration_results.tex'))
}

{
  
  Covariance <- function( h, theta, delta_taper, ... ) {
    return(cov.wend1_1(h, theta, ...) * cov.wend2_1(h, c(delta_taper, 1, 0), ...))
  }
  
  delta_taper <- 6
  tikz(paste0(CommonPath, '/finalplots/covariograms_results.tex'), standAlone = T, width = 7.5, height = 3.75)
  par(pty = 'm', oma = c(0, 2.5, 0, 0), mar = c(2, 0, 1.5, 1), mfrow = c(1, 1))
  plot(x = seq(0, 20, length.out = 2000), xaxs = 'i', axes = FALSE, ann = FALSE,
       y = cov.wend1_1(seq(0, 20, length.out = 2000), theta = output_taper$par[1:3]), 
       type = 'l', ylim = c(0, 0.2) , xlim = c(0, 20), col=colors_to_use[2], lty = 1, lwd = 7)
  lines(x = seq(0, 20, length.out = 4000),
        y = cov.wend2_1(seq(0, 20, length.out = 4000), theta = output_geofit$param[c(2, 3, 1)]), col = colors_to_use[3], lwd = 7)
  
  lines(x = seq(0, 20, length.out = 4000),
        y = cov.wend2_1(seq(0, 20, length.out = 4000), theta = c(delta_taper, output_directmiss$par[1:2])), col = colors_to_use[4], lwd = 7 )
  
  lines(x = seq(0, 20, length.out = 4000),
        y = Covariance(seq(0, 20, length.out = 4000), theta = output_taper$par[1:3], delta_taper = delta_taper), type = 'l', 
        col = colors_to_use[2], lty = 2, lwd = 7)
  
  legend('topright', bty = 'n', legend = c('CL', 'TA estimated cov.', 'TA tapered cov.', 'DM'),
         col = c(colors_to_use[3], colors_to_use[2], colors_to_use[2], colors_to_use[4]), lty = c(1, 1, 2, 1) , lwd = c(7, 7, 7, 7), cex = 1.75)
  box()
  axis(1, cex.axis = 2)
  axis(2, cex.axis = 2)
  dev.off()
  tinytex::pdflatex(paste0(CommonPath, '/finalplots/covariograms_results.tex'))
  system(paste0('rm ~', '/finalplots/covariograms_results.tex'))
}
