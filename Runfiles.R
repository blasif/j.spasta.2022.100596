
# Install packages -----

rm(list = ls())
system('Rscript docker_run.R')

# Simulation study -----

# Create Fields

rm(list = ls())
source('00_system.R')
system('Rscript 02_create_fields.R')

# Estimation Part

list_fields <- list.files(paste0(CommonPath, '/fields'))
list_fields <- list_fields[grep('field.', list_fields)]

system('Rscript --version') # !

for(i in 1:length(list_fields)){
  
  print(i)
  
  # Classical ML
  
  system(paste0(RScript, ' 03a_estim_matern.R', ' ', list_fields[i]))
  
  # Taper approach
  
  system(paste0(RScript, ' 03b_estim_taper.R', ' ', list_fields[i]))
  
  # Direct misspecification
  
  system(paste0(RScript, ' 03c_estim_directmiss.R', ' ', list_fields[i]))
  
  # Composite likelihood
  
  system(paste0(RScript, ' 03d_estim_composite.R', ' ', list_fields[i]))
  
}

# Simulation figures

rm(list = ls())
source('00_system.R')

system(paste0(RScript, ' 04_figures.R'))

# Simulation figures not in the article

system(paste0(RScript, ' 04_figures_extra.R'))

# Application -----

system(paste0(RScript, ' ~/application/extractNetCDF.R'))

system(paste0(RScript, ' ~/application/illustration.R'))

# Application figures

system(paste0(RScript, ' ~/application/illu_figures.R'))
