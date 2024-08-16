rm(list = ls())

source('00_system.R')

# Create the baseline simulation profiles

if(create_new <- T){
  
  df_settings <- expand.grid('parsetting' = c('low', 'high'), 'gridded' = c(T, F),
                             'nsize' = c(784L, 3025L), 'gaussian' = c(T, F), stringsAsFactors = F)
  settings_list <- list()
  for(i in 1:dim(df_settings)[1]){
    settings_list[[i]] <- new(Class = 'fieldspec', parsetting = df_settings[i, 1],
                              gridded = df_settings[i, 2], nsize = df_settings[i, 3],
                              gaussian = df_settings[i, 4], specsnr = i)
  }
  
  save(settings_list, file = paste0(CommonPath, '/fields/settings.RData'))
} else(load(paste0(CommonPath, '/fields/settings.RData')))

# Validate simulation profiles (this will simulate the fields and fill all the empty slots for the class fieldspec)

for(zz in 1:length(settings_list)){
  
  this_settings <- settings_list[[zz]]
  
  { cl <- makeCluster(ncores-2)
    setDefaultCluster(cl = cl)
    clusterExport(cl, c('lib.loc'))
    clusterEvalQ(cl, .libPaths(new = lib.loc))
    clusterEvalQ(cl, library('spam', lib.loc = lib.loc))
    clusterEvalQ(cl, library('spam64', lib.loc = lib.loc))
    clusterEvalQ(cl, library('spatstat', lib.loc = lib.loc))
    clusterEvalQ(cl, library('GeoModels', lib.loc = lib.loc))
    clusterExport(cl, c( 'zz', 'this_settings' , 'CommonPath', 'seed_master_number', 'cov.mat2'))
    
    parSapply(cl, 1:this_settings@fielditer_n, FUN = function(iter) {
      
      source('01_classes.R')
      
      print(paste0('setting: ', zz, '.iter: ', iter))
      
      this_field <- validate(theObject = this_settings, i_ref = iter)
      
      save(this_field, file = paste0(CommonPath, '/fields/field.', zz, '.iter', iter, '.RData'))
    })
  }; stopCluster(cl)
} 
