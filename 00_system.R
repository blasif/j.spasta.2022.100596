
# Main variables & more

CommonPath <- getwd( )                                     # main path project
RScript <- 'Rscript'                                      # Rscript command to run scripts
lib.loc <- paste0(CommonPath, '/lib.loc')                  # own library
.libPaths(new = lib.loc)                                  # specifying new library           
seed_master_number <- 1778                                # Master seed number. This influences all other seed in the code

if (version$major != '4'){
  stop('use R-4.x.y')
} 

# load packages

{
  suppressMessages(require(spatstat, lib.loc = lib.loc, quietly = TRUE, warn.conflicts = FALSE))
  suppressMessages(require(spam, lib.loc = lib.loc, quietly = TRUE, warn.conflicts = FALSE))
  require(spam64, lib.loc = lib.loc, quietly = TRUE, warn.conflicts = FALSE)
  suppressMessages(require(devtools, lib.loc = lib.loc, quietly = TRUE, warn.conflicts = FALSE))
  require(optimParallel, lib.loc = lib.loc, quietly = TRUE, warn.conflicts = FALSE)
  suppressMessages(require(GeoModels, lib.loc = lib.loc, quietly = TRUE, warn.conflicts = FALSE))
  suppressMessages(require(vioplot, lib.loc = lib.loc , quietly = TRUE, warn.conflicts = FALSE))
  require(latex2exp, lib.loc = lib.loc, quietly = TRUE)
  require(tinytex, lib.loc = lib.loc , quietly = TRUE)
  options('tikzLatex' = '/home/rstudio/.TinyTeX/bin/x86_64-linux/pdflatex')
  require(tikzDevice, lib.loc = lib.loc, quietly = TRUE)
  require(numDeriv, lib.loc = lib.loc, quietly = TRUE)
  suppressMessages(require(mgcv, lib.loc = lib.loc, quietly = TRUE))
}

# call other required R scripts

source('01_classes.R', verbose = FALSE)
source('99_helper_functions.R', verbose = FALSE)

# detect ncores

ncores <- parallel::detectCores()
