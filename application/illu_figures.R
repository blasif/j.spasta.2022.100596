rm(list = ls())

source('00_system.R')

colors_to_use <- palette.colors(palette = 'Okabe-Ito')

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
