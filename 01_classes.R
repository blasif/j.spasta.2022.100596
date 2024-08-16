
# Fields class -------

setClass('fieldspec',
         # define slots and respective types
         slots = list(
           
           # Admin -------------------- #
           specsnr      = 'integer' ,       # specs row of /fields/settings.RData
           parsetting   = 'character' ,     # high/low
           gaussian     = 'logical' ,       # is Gaussian? T/F
           extra_par    = 'vector' ,        # non-gaussian field extra-parameters
           fielditer_n  = 'integer' ,       # number of iterations for the same field. Preset to 200
           fielditer    = 'integer' ,       # iteration number
           seed         = 'integer' ,       # seed used to generate this field (master_seed_number * fielditer)

           # Field parameters --------- #
           nsize     = 'integer' ,          # size of the field. Filled later on based on parsetting
           gridded   = 'logical' ,          # if TRUE = gridded , FALSE = ssh
           min_dist  = 'numeric' ,          # if gridded = FALSE, min distance between locations
           domain    = 'vector' ,           # 2x2 matrix with boundaries for the domain
           truetheta = 'vector' ,           # true 4 length vector thetas (range, psill, smoothness, nugget)
           truebetas = 'vector' ,           # true 3 length vector of betas (beta0, beta_x, beta_y)
           
           # Store info --------------- #
           locs            = 'matrix' ,      # expand.grid output
           fulldata        = 'vector'       # complete matern simulation observations
           ),
         
         prototype = list(
           
           # Admin -------------------- #
           specsnr      = NA_integer_ ,
           parsetting   = NA_character_ , 
           gaussian     = TRUE , 
           extra_par    = vector() , 
           fielditer_n  = 200L , 
           fielditer    = NA_integer_ , 
           seed         = NA_integer_ , 

           # Field parameters --------- #
           nsize     = NA_integer_ , 
           gridded   = T , 
           min_dist  = 0.01 , 
           domain    = matrix(c(0,1, 0,1), nrow = 2, ncol = 2, byrow = TRUE) ,
           truetheta = vector() , 
           truebetas = c(1, 0.1, 0.2) , 
           
           # store info --------------- #
           locs            = matrix() , 
           fulldata        = vector()  
         )
)

# Define validation function for the Field class ----

setGeneric(name = 'validate',
           def = function(theObject, i_ref) {
             standardGeneric('validate') })

setMethod(f = 'validate',
          signature = 'fieldspec',
          definition = function(theObject, i_ref) {
            
            # update seed and rep
            
            theObject@fielditer <- as.integer(i_ref)
            theObject@seed <- as.integer(i_ref * seed_master_number)
            
            set.seed(seed = theObject@seed) # reproducibility
            
            # par_setting
            
            if(theObject@parsetting == 'low'){
              theObject@truetheta <-  c(0.05, 1, 0.5, 0.1)
            } else { 
              if(theObject@parsetting == 'high'){
                theObject@truetheta <-  c(0.05, 1, 2.3, 0)
                } else stop('parsetting argument not correct')
              }
            
            # sampling_design design
            
            if(theObject@gridded){ # Gridded
              
              x_seq  <- seq(theObject@domain[1,1], theObject@domain[1,2], length.out = sqrt(theObject@nsize))
              y_seq  <- seq(theObject@domain[2,1], theObject@domain[2,2], length.out = sqrt(theObject@nsize))
              theObject@locs <- as.matrix(expand.grid(x = x_seq, y = y_seq))

            } else{ # SSI
              
              aux_domain <- owin(c(theObject@domain[1,1], theObject@domain[1,2]), c(theObject@domain[2,1], theObject@domain[2,2]))
              aux_rSSI <- rSSI(r = theObject@min_dist, n = theObject@nsize, win = aux_domain)
              locs <- cbind('x_seq' = aux_rSSI[[3]],'y_seq' = aux_rSSI[[4]])
              theObject@locs <- locs[order( locs[,1]), ]
              
            }
            
            # simulate field
            
            if (theObject@gaussian){ # then gaussian
              
              # complete distance matrix
              distmat <- as.matrix(dist(theObject@locs))
              Sigma <- cov.mat2(h = distmat, theta = theObject@truetheta)
              
              rm(distmat)
              # gc(verbose = FALSE, reset = T)

              iidsample <- rnorm(theObject@nsize)
              cholS <- chol(Sigma)
              
              partial <- theObject@truebetas[1] + theObject@truebetas[2] * cos(theObject@locs[,1]) + 
                theObject@truebetas[3]*cos(theObject@locs[,2])
              
              theObject@fulldata <- partial + as.vector(iidsample %*% cholS)
              
            } else {
              
              theObject@extra_par <- c(0.5, 0.8) # skew, tail (kurtosis)
              
              # Simulation of the spatial SAS RF:
              tmp_data <- GeoSim(coordx = theObject@locs, corrmodel = 'Matern', model = 'SinhAsinh', 
                             param = list(skew = theObject@extra_par[1], tail = theObject@extra_par[2], 
                                         smooth = theObject@truetheta[3], mean = 0, 
                                         sill = theObject@truetheta[2], scale = theObject@truetheta[1],
                                         nugget = theObject@truetheta[4]))$data
              
              # Outputs
              
              partial <- theObject@truebetas[1] + theObject@truebetas[2] * cos(theObject@locs[,1]) + 
                theObject@truebetas[3]*cos(theObject@locs[,2])
              
              theObject@fulldata <- partial + tmp_data
              
            }
            
            return(theObject)
            
          }
)

# Estimation classes --------

# 03a full matern class -----

setClass('matern_estimate',
         # define slots and respective types
         slots = list(
           
           # admin
           specsnr = 'integer',
           fielditer = 'integer',
           
           # subsample
           time = 'numeric',
           optimout = 'ANY',
           betas_hat = 'vector',
           memory_used = 'numeric',
           loadstring = 'character'),
         
         prototype = list(
           
           # admin
           specsnr = NA_integer_,
           fielditer = NA_integer_,
           
           # subsample
           time = NA_real_,
           optimout = NULL,
           betas_hat = vector(),
           memory_used = NA_real_,
           loadstring = NA_character_)
)

# 03b tapering class -----

setClass('taper_estimate',
         # define slots and respective types
         slots = list(
           
           # admin
           specsnr = 'integer',
           fielditer = 'integer',
           delta = 'numeric',
           
           # subsample
           time = 'numeric',
           optimout = 'ANY',
           betas_hat = 'vector',
           memory_used = 'numeric',
           loadstring = 'character'),
         
         prototype = list(
           
           # admin
           specsnr = NA_integer_,
           fielditer = NA_integer_,
           delta = NA_real_,
           
           # subsample
           time = NA_real_,
           optimout = NULL,
           betas_hat = vector(),
           memory_used = NA_real_,
           loadstring = NA_character_)
)

# 03c directmiss class -----

setClass('directmiss_estimate',
         # define slots and respective types
         slots = list(
           
           # admin
           specsnr = 'integer',
           fielditer = 'integer',
           delta = 'numeric',
           
           # subsample
           time = 'numeric',
           optimout = 'ANY',
           betas_hat = 'vector',
           memory_used = 'numeric',
           loadstring = 'character'),
         
         prototype = list(
           
           # admin
           specsnr = NA_integer_,
           fielditer = NA_integer_,
           delta = NA_real_,
           
           # subsample
           time = NA_real_,
           optimout = NULL,
           betas_hat = vector(),
           memory_used = NA_real_,
           loadstring = NA_character_)
)

# 03d composite class -----

setClass('composite_estimate',
         # define slots and respective types
         slots = list(
           
           # admin
           specsnr = 'integer',
           fielditer = 'integer',
           delta = 'numeric',
           
           # subsample
           time = 'numeric',
           optimout = 'ANY',
           betas_hat = 'vector',
           memory_used = 'numeric',
           loadstring = 'character'),
         
         prototype = list(
           
           # admin
           specsnr = NA_integer_,
           fielditer = NA_integer_,
           delta = NA_real_,
           
           # subsample
           time = NA_real_,
           optimout = NULL,
           betas_hat = vector(),
           memory_used = NA_real_,
           loadstring = NA_character_)
)
