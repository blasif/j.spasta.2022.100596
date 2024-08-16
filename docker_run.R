CommonPath <- getwd()                                     # main path project
RScript <- 'Rscript'                                      # Rscript command to run scripts
lib.loc <- paste0(CommonPath, '/lib.loc')                  # own library
.libPaths(new = lib.loc)                                  # specifying new library           

install_package_source <- function(packagename, version, lib = lib.loc){
  
  # search in archive
  
  pkg_name <- paste0('https://cran.r-project.org/src/contrib/Archive/',
                     packagename, '/', packagename, '_', version, '.tar.gz')
  
  check <- suppressWarnings(try(open.connection(url(pkg_name), open = 'rt'), silent = T)[1])
  
  # install from archive
  if(ifelse(is.null(check), T, F)){
    
    install.packages(pkg_name, repos = NULL, type = 'source', dependencies = NA, lib = lib.loc)
    
    close(url(pkg_name))
  } else {
    
    # install current version
    install.packages(paste0('https://cran.r-project.org/src/contrib/', packagename, '_', version, '.tar.gz'), repos = NULL, type = 'source', dependencies = NA, lib = lib.loc)
    close(url(paste0('https://cran.r-project.org/src/contrib/', packagename, '_', version, '.tar.gz')))
  }
  
  ifelse(length(installed.packages(lib.loc = lib.loc)[which(installed.packages(lib.loc = lib.loc)[, 1] == packagename)]) == 0, return(F), return(T))
  
}

CommonPath <- getwd()
lib.loc <- paste0(CommonPath, '/lib.loc')
.libPaths(new = lib.loc)

pkgs_to_install <- read.csv(file = 'pkgs.csv')[, c(2, 3)]

pkgs_to_install <- cbind.data.frame(pkgs_to_install, 'status' = logical(length = dim(pkgs_to_install)[1]))

Geomodels_info <- pkgs_to_install[26, ]

pkgs_to_install <- pkgs_to_install[c(-26, -157), ]

for(i in 1:dim(pkgs_to_install)[1]){
  
  pkgs_to_install[i, 3] <- install_package_source(packagename = pkgs_to_install[i, 1], 
                                                  version = pkgs_to_install[i, 2])
  
}

remotes::install_github('vmoprojs/GeoModels', ref = 'v1.1',dependencies = FALSE, lib = lib.loc)

ifelse(tinytex::is_tinytex(), print('tinytex installed'), tinytex::install_tinytex())

tinytex::tlmgr_install('pgf')
tinytex::tlmgr_install('grfext')
tinytex::tlmgr_install('preview')