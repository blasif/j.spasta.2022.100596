Data via https://esgf-data.dkrz.de/search/cmip6-dkrz/ based on filter  
CMIP6 | CESM2 | tas | abrupt-4xCO2 | r1i1p1f1 | Amon  
select first file. Direct access  
http://esgf-data.ucar.edu/thredds/fileServer/esg_dataroot/CMIP6/CMIP/NCAR/CESM2/abrupt-4xCO2/r1i1p1f1/Amon/tas/gn/v20190927/tas_Amon_CESM2_abrupt-4xCO2_r1i1p1f1_gn_000101-015012.nc

The file is copied to `$(pwd)/application`. It contains monthly means (mid month) for 150 years with a four fold CO2 emission.
Have a glance with 
```
ncview $(pwd)/tas_Amon_* &
```
then clicking `tas` and then fast forward button.


Reading the data is done with the R script in this directory.