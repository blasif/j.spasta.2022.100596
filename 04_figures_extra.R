rm(list = ls())

source('00_system.R')

colors_to_use <- palette.colors(palette = 'Okabe-Ito')
do_analysis <- F

list_estimates <- list.files(paste0(CommonPath, '/estimates'))
list_estimates <- list_estimates[grep('estimates.', list_estimates)]

to_fill <- matrix(NA, ncol = 4, nrow = 200)
colnames(to_fill) <- c('Matern', 'Tapering', 'Directmiss', 'Composite')

# SAS parameters -----

result_kurt <- lapply(1:16, function(x) to_fill)
result_skew <- lapply(1:16, function(x) to_fill)

bias_analysis_kurt <- data.frame('bias' = NA, 'gridded' = NA, 'gauss' = NA, 'size' = NA, 'thetas' = NA, 'estim'= NA)
bias_analysis_skew <- data.frame('bias' = NA, 'gridded' = NA, 'gauss' = NA, 'size' = NA, 'thetas' = NA, 'estim'= NA)

for(i in 1:length(list_estimates)){
  
  load(paste0(CommonPath, '/estimates/', list_estimates[i]))
  load(paste0(CommonPath, '/fields/', this_estimates@loadstring))
  
  if(this_field@gaussian) next
  
  if(is(this_estimates, 'composite_estimate')){
    
    z <- 4
    
    result_kurt[[this_estimates@specsnr]][this_estimates@fielditer, z] <- unlist(this_estimates@optimout$par[8])
    result_skew[[this_estimates@specsnr]][this_estimates@fielditer, z] <- unlist(this_estimates@optimout$par[7])
    
  } else {
    
    if(is(this_estimates, 'directmiss_estimate')){
      
      z <- 3
      result_kurt[[this_estimates@specsnr]][this_estimates@fielditer, z] <- this_estimates@optimout$par[4]
      result_skew[[this_estimates@specsnr]][this_estimates@fielditer, z] <- this_estimates@optimout$par[3]
      
    } else{
      
      if(is(this_estimates, 'taper_estimate')){z <- 2}
      if(is(this_estimates, 'matern_estimate')){z <- 1}
      
      result_kurt[[this_estimates@specsnr]][this_estimates@fielditer, z] <- this_estimates@optimout$par[5]
      result_skew[[this_estimates@specsnr]][this_estimates@fielditer, z] <- this_estimates@optimout$par[4]
      
    }
    
  }
  
  bias_analysis_kurt[i, ] <- c(result_kurt[[this_estimates@specsnr]][this_estimates@fielditer, z], this_field@gridded, this_field@gaussian, this_field@nsize, this_field@parsetting, z)
  bias_analysis_skew[i, ] <- c(result_skew[[this_estimates@specsnr]][this_estimates@fielditer, z], this_field@gridded, this_field@gaussian, this_field@nsize, this_field@parsetting, z)
  
}

# Betas ----

result_beta_1 <- lapply(1:16, function(x) to_fill)
result_beta_2 <- lapply(1:16, function(x) to_fill)

bias_analysis_beta_1 <- data.frame('bias' = NA_real_, 'gridded' = NA, 'gauss' = NA, 'size' = NA_integer_, 'thetas' = NA, 'estim'= NA)
bias_analysis_beta_2 <- data.frame('bias' = NA_real_, 'gridded' = NA, 'gauss' = NA, 'size' = NA_integer_, 'thetas' = NA, 'estim'= NA)

for(i in 1:length(list_estimates)){
  load(paste0(CommonPath, '/estimates/', list_estimates[i]))
  load(paste0(CommonPath, '/fields/', this_estimates@loadstring))
  
  if(is(this_estimates, 'composite_estimate')){z <- 4}
  if(is(this_estimates, 'directmiss_estimate')){z <- 3}
  if(is(this_estimates, 'taper_estimate')){z <- 2}
  if(is(this_estimates, 'matern_estimate')){z <- 1}
  
  result_beta_1[[this_estimates@specsnr]][this_estimates@fielditer, z] <- this_estimates@betas_hat[2]
  result_beta_2[[this_estimates@specsnr]][this_estimates@fielditer, z] <- this_estimates@betas_hat[3]
  
  if(i == 1){
    
    bias_analysis_beta_1 <- data.frame('bias' = (as.numeric(result_beta_1[[this_estimates@specsnr]][this_estimates@fielditer, z]-this_field@truebetas[2])/this_field@truebetas[2]),
                                       'gridded' = this_field@gridded, 'gauss' = this_field@gaussian, 'size' = factor(as.character(this_field@nsize), levels = c('784', '3025')),
                                       'thetas' = factor(this_field@parsetting, levels = c('low', 'high')), 'estim' = factor(as.character(z), levels = c('1', '2', '3', '4')))
    
    bias_analysis_beta_2 <- data.frame('bias' = (as.numeric(result_beta_2[[this_estimates@specsnr]][this_estimates@fielditer, z]-this_field@truebetas[3])/this_field@truebetas[3]),
                                       'gridded' = this_field@gridded, 'gauss' = this_field@gaussian, 'size' = factor(as.character(this_field@nsize), levels = c('784', '3025')),
                                       'thetas' = factor(this_field@parsetting, levels = c('low', 'high')), 'estim' = factor(as.character(z), levels = c('1', '2', '3', '4')))
    
  } else {
    
    bias_analysis_beta_1 <- rbind.data.frame(bias_analysis_beta_1, setNames(
      data.frame((as.numeric(result_beta_1[[this_estimates@specsnr]][this_estimates@fielditer, z]-this_field@truebetas[2])/this_field@truebetas[2]),
                 this_field@gridded,
                 this_field@gaussian,
                 this_field@nsize,
                 this_field@parsetting,
                 z), names(bias_analysis_beta_1)))
    
    bias_analysis_beta_2 <- rbind.data.frame(bias_analysis_beta_2, setNames(
      data.frame((as.numeric(result_beta_2[[this_estimates@specsnr]][this_estimates@fielditer, z]-this_field@truebetas[3])/this_field@truebetas[3]),
                 this_field@gridded,
                 this_field@gaussian,
                 this_field@nsize,
                 this_field@parsetting,
                 z), names(bias_analysis_beta_2)))
    
  }
}

# Beta 1 ----

if(do_analysis){
  
  plot(density(bias_analysis_beta_1$bias[complete.cases(bias_analysis_beta_1$bias)]))
  model_beta_1 <- lm(bias ~ . + estim*thetas + estim*gridded + estim*gauss + estim*size, data = bias_analysis_beta_1)
  anova(model_beta_1)
  plot(model_beta_1)
}

{
  
  mean_values_beta_1 <- aggregate(bias~estim, data = bias_analysis_beta_1, FUN = mean)$'bias'
  mean_values_beta_1_b <- aggregate(bias~estim + thetas, data = bias_analysis_beta_1, FUN = mean)$'bias'
  boxplot(data = bias_analysis_beta_1, bias ~  gauss + size + gridded +  thetas + estim)
  
  text_names <- c('beta1_a.tex', 'beta1_b.tex', 'beta1_c.tex', 'beta1_d.tex')
  
  lower_bound <- -35
  
  for(i in 1:4){
    
    tikz(paste0(CommonPath, '/extraplots/', text_names[i]), standAlone = T, width = 4, height = 6)
    
    par(pty = 'm', oma = c(3, 3, 0, 0), mar = c(2, 0, 1.5, 0.5), mfrow = c(1, 1))
    plot(0, type = 'n', axes = FALSE, ann = FALSE, ylim= c(lower_bound, 35), xlim = c(0.5, 4.5))
    axis(side = 2, cex.axis = 2, tck = 1, lty = 2, col = 'grey70', labels = ifelse(i == 1, T, F))
    box()
    
    boxplot(data = bias_analysis_beta_1[which(bias_analysis_beta_1$estim == i), ], bias ~  gauss + thetas, 
            axes = FALSE, cex.box = 2, boxlwd = 4, outwex = 3, boxwex = 0.75, add = T, main = '')
    
    if(i == 1) {
      segments(0.5, mean_values_beta_1[1], 4.5, mean_values_beta_1[1], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_beta_1_b[1], 2.5, mean_values_beta_1_b[1], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_beta_1_b[5], 4.5, mean_values_beta_1_b[5], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    if(i == 2){
      segments(0.5, mean_values_beta_1[2], 4.5, mean_values_beta_1[2], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_beta_1_b[2], 2.5, mean_values_beta_1_b[2], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_beta_1_b[6], 4.5, mean_values_beta_1_b[6], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    if(i == 3){
      segments(0.5, mean_values_beta_1[2], 4.5, mean_values_beta_1[2], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_beta_1_b[3], 2.5, mean_values_beta_1_b[3], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_beta_1_b[7], 4.5, mean_values_beta_1_b[7], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    if(i == 4){
      segments(0.5, mean_values_beta_1[3], 4.5, mean_values_beta_1[3], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_beta_1_b[4], 2.5, mean_values_beta_1_b[4], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_beta_1_b[8], 4.5, mean_values_beta_1_b[8], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    
    # first line labels
    axis(side = 1, at = c(1, 2, 3, 4), cex.axis = 2, labels = c('SAS', 'Gau', 'SAS', 'Gau'))
    
    # second line labels
    
    axis(side = 1, at = c(1, 1.5, 2), cex.axis = 2, labels = c(NA, NA, NA), outer = T, tck = 0, pos = lower_bound - 10)
    
    axis(side = 1, at = c(1.5), cex.axis = 2, labels = c('Low'), outer = T, tck = F, pos = lower_bound - 10)
    
    axis(side = 1, at = c(1.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos = lower_bound - 10)
    axis(side = 1, at = c(1, 2), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound - 10)
    
    axis(side = 1, at = c(3, 3.5, 4), cex.axis = 2, labels = c(NA, NA, NA), outer = T, tck = 0, pos = lower_bound - 10)
    axis(side = 1, at = c(3.5), cex.axis = 2, labels = c('High'), outer = T, tck = F, pos = lower_bound - 10)
    axis(side = 1, at = c(3.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos = lower_bound - 10)
    axis(side = 1, at = c(3,4), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound - 10)
    
    
    dev.off()
    
    tinytex::pdflatex(paste0(CommonPath, '/extraplots/', text_names[i]))
    
    system(paste0('rm ~', '/extraplots/', text_names[i]))
    
  }
}

# Beta 2 ----

if(do_analysis){
  
  plot(density(bias_analysis_beta_2$bias[complete.cases(bias_analysis_beta_2$bias)]))
  model_beta_2 <- lm(bias ~ . + estim*thetas + estim*gridded + estim*gauss + estim*size, data = bias_analysis_beta_2)
  anova(model_beta_2)
  plot(model_beta_2)
}

{
  
  mean_values_beta_2 <- aggregate(bias~estim, data = bias_analysis_beta_2, FUN = mean)$'bias'
  mean_values_beta_2_b <- aggregate(bias~estim + thetas, data = bias_analysis_beta_2, FUN = mean)$'bias'
  boxplot(data = bias_analysis_beta_2, bias ~  gauss + size + gridded +  thetas + estim)
  
  text_names <- c('beta2_a.tex', 'beta2_b.tex', 'beta2_c.tex', 'beta2_d.tex')
  
  lower_bound <- -35
  
  for(i in 1:4){
    
    tikz(paste0(CommonPath, '/extraplots/', text_names[i]), standAlone = T, width = 4, height = 6)
    
    par(pty = 'm', oma = c(3, 3, 0, 0), mar = c(2, 0, 1.5, 0.5), mfrow = c(1, 1))
    plot(0, type = 'n', axes = FALSE, ann = FALSE, ylim= c(lower_bound, 35), xlim = c(0.5, 4.5))
    axis(side = 2, cex.axis = 2, tck = 1, lty = 2, col = 'grey70', labels = ifelse(i == 1, T, F))
    box()
    
    boxplot(data = bias_analysis_beta_2[which(bias_analysis_beta_2$estim == i), ], bias ~  gauss + thetas, 
            axes = FALSE, cex.box = 2, boxlwd = 4, outwex = 3, boxwex = 0.75, add = T, main = '')
    
    if(i == 1) {
      segments(0.5, mean_values_beta_2[1], 4.5, mean_values_beta_2[1], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_beta_2_b[1], 2.5, mean_values_beta_2_b[1], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_beta_2_b[5], 4.5, mean_values_beta_2_b[5], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    if(i == 2){
      segments(0.5, mean_values_beta_2[2], 4.5, mean_values_beta_2[2], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_beta_2_b[2], 2.5, mean_values_beta_2_b[2], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_beta_2_b[6], 4.5, mean_values_beta_2_b[6], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    if(i == 3){
      segments(0.5, mean_values_beta_2[2], 4.5, mean_values_beta_2[2], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_beta_2_b[3], 2.5, mean_values_beta_2_b[3], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_beta_2_b[7], 4.5, mean_values_beta_2_b[7], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    if(i == 4){
      segments(0.5, mean_values_beta_2[3], 4.5, mean_values_beta_2[3], col = colors_to_use[2], lty = 1, lwd = 6)
      segments(0.5, mean_values_beta_2_b[4], 2.5, mean_values_beta_2_b[4], col = colors_to_use[6], lty = 2, lwd = 4)
      segments(2.5, mean_values_beta_2_b[8], 4.5, mean_values_beta_2_b[8], col = colors_to_use[6], lty = 2, lwd = 4)
    }
    
    # first line labels
    axis(side = 1, at = c(1, 2, 3, 4), cex.axis = 2, labels = c('SAS', 'Gau', 'SAS', 'Gau'))
    
    # second line labels
    
    axis(side = 1, at = c(1, 1.5, 2), cex.axis = 2, labels = c(NA, NA, NA), outer = T, tck = 0, pos = lower_bound - 40)
    
    axis(side = 1, at = c(1.5), cex.axis = 2, labels = c('Low'), outer = T, tck = F, pos = lower_bound - 10)
    
    axis(side = 1, at = c(1.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos = lower_bound - 10)
    axis(side = 1, at = c(1, 2), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound - 10)
    
    axis(side = 1, at = c(3, 3.5, 4), cex.axis = 2, labels = c(NA, NA, NA), outer = T, tck = 0, pos = lower_bound - 40)
    axis(side = 1, at = c(3.5), cex.axis = 2, labels = c('High'), outer = T, tck = F, pos = lower_bound - 10)
    axis(side = 1, at = c(3.5), cex.axis = 2, labels = c(NA), outer = T, tck = -0.01, pos = lower_bound - 10)
    axis(side = 1, at = c(3,4), cex.axis = 2, labels = c(NA), outer = T, tck = 0.01, pos = lower_bound - 10)
    
    
    dev.off()
    
    tinytex::pdflatex(paste0(CommonPath, '/extraplots/', text_names[i]))
    
    system(paste0('rm ~', '/extraplots/', text_names[i]))
    
  }
}

# Kurtosis ----

bias_analysis_kurt <- bias_analysis_kurt[!is.na(bias_analysis_kurt$gauss), ]
bias_analysis_kurt$bias  <- as.numeric(bias_analysis_kurt$bias)
bias_analysis_kurt$estim <- as.factor(bias_analysis_kurt$estim)

if(do_analysis){
  
  plot(density(bias_analysis_kurt$bias))
  plot(density(log(bias_analysis_kurt$bias)))
  model_kurt <- lm(log(bias) ~ gridded + size + thetas + estim , data = bias_analysis_kurt)
  anova(model_kurt)
  plot(model_kurt)
  
}

pdf(paste0(CommonPath, '/extraplots/', 'kurt.pdf'), width = 7.5, height = 8)
par(pty = 'm', oma = c(0, 0, 0, 0), mar = c(2, 2.5, 1.5, 2.5))
boxplot(data = bias_analysis_kurt, bias ~ estim, cex.axis = 1.5, cex.lab = 2, axes = FALSE, 
        cex.box = 2, boxlwd = 3, outwex = 1, boxwex = 0.75)
box()
axis(side = 1, at = c(1, 2, 3, 4), cex.axis = 1, labels = c('ML', 'TA', 'DM', 'CL'))
axis(side = 2)
abline(h = 0.8, lty = 2, col = 'blue')
grid(NA, 7)

dev.off()

# Skewness ----

bias_analysis_skew <- bias_analysis_skew[!is.na(bias_analysis_skew$gauss), ]
bias_analysis_skew$bias  <- as.numeric(bias_analysis_skew$bias)
bias_analysis_skew$estim <- as.factor(bias_analysis_skew$estim)

if(do_analysis){
  
  plot(density(bias_analysis_skew$bias))
  plot(density(log(bias_analysis_skew$bias)))
  model_skew <- lm(bias ~ gridded + size + thetas + estim , data = bias_analysis_skew)
  anova(model_skew)
  model_skew
  plot(model_skew)
  
}

pdf(paste0(CommonPath, '/extraplots/', 'skew.pdf'), width = 7.5, height = 8)
par(pty = 'm', oma = c(0, 0, 0, 0), mar = c(2, 2.5, 1.5, 2.5))
boxplot(data = bias_analysis_skew, bias ~ estim, cex.axis = 1.5, cex.lab = 2, axes = FALSE, 
        cex.box = 2, boxlwd = 3, outwex = 1, boxwex = 0.75)
box()
axis(side = 1, at = c(1, 2, 3, 4), cex.axis = 1, labels = c('ML', 'TA', 'DM', 'CL'))
axis(side = 2)
grid(NA, 7)
abline(h = 0.5, lty = 2, col = 'blue')
dev.off()
