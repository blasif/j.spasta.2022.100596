
rm(list = ls())
source('00_system.R')

ifelse(length(which(installed.packages()[,1] == 'ncdf4')) == 0 , install.packages('ncdf4',lib = lib.loc), print('ncdf4 installed'))

require(ncdf4)  # installation might be necessary

path <- paste0(CommonPath,'/application/')
file <- 'tas_Amon_CESM2_abrupt-4xCO2_r1i1p1f1_gn_000101-015012.nc'

nc <- nc_open(paste0(path, file))

# just for illustration, can be seen from the `ncview` tool
nc$nvars

long <- ncvar_get(nc, 'lon') 
lati <- ncvar_get(nc, 'lat')
tas <- ncvar_get(nc, 'tas')

# consider JJA, years 131:150 minus 1:20
jja2 <- 6:8 + rep((130:149)*12, each = 3)
jja1 <- 6:8 + rep((0:19)*12, each = 3)

# extract data out of the three-dim array and average:

tas1 <- apply( tas[, , jja1], c(1, 2), mean)
tas2 <- apply( tas[, , jja2], c(1, 2), mean)
delta <- tas2 - tas1 

require(fields,lib.loc = lib.loc)
require(maps, lib.loc = lib.loc)
image.plot(long, lati, tas1)
map(add = TRUE, wrap = c(0, 360))

image.plot(long, lati, delta)
map(add = TRUE, wrap = c(0, 360))

# We have to possibly eliminate top and bottom row due to numerical instabilities in the great circle distance.
save(tas1, tas2, lati, long, file = paste0(path,'ClimChange.RData'))

nc_close(nc)
rm(tas)
gc()
