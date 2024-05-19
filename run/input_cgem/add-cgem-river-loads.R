#-- add-tp-river-loads.R
#L3, 2023
# Starts with 2013_leem_fine_river_data.nc

#To read netCDF model file
library(ncdf4)

setwd("~/FVCOM_CGEM")
source("NUTRIENT_CGEM.R")

# Location of Original 'just hydro' river forcing file
#hydrofile <- "/Users/lllowe/LakeErie/TPmodel/simulations/input/full_year/2013_leem_fine_river_data_original_hydro.nc"
#Or used cropped version:
hydrofile <- "/Users/lllowe/LakeErie/TPmodel/simulations/input/2013_leem_fine_river_data.nc"
# New river file with added TP loads
tpfile <- "/Users/lllowe/FVCOM_CGEM/2013_leem_fine_river_data.nc"

#Check file paths
if (!file.exists(hydrofile)) stop("Original river loads file not found, exiting")

# Copy the original netCDF file, since we will append to the file
file.copy(hydrofile,tpfile,overwrite=TRUE)

#read the file
ncdata <- nc_open(tpfile,write=TRUE,readunlim=FALSE)

# Get Data, Times (MJD, days since 1858-11-17 00:00:00)
time <- ncvar_get(ncdata, "Times", start=c(1,1), count=c(-1,-1))
time <- as.POSIXct(time, format="%Y-%m-%dT%H:%M:%S", tz="GMT")
nt <- length(time)
rivnames <- ncvar_get(ncdata, "river_names", start=c(1,1), count=c(-1,-1))
nr <- length(rivnames)

#Start making the new array to put in netCDF file
# Define netCDF dimensions
dimrivers <- ncdata$dim[['rivers']]
dimtime <- ncdata$dim[['time']]

#Used for missing value
mv2 <- 1.e30

# Define the new variables
for(i in 1:length(cgem_vars)){
  value <- ncvar_def( name=cgem_vars[i], units=cgem_units[i], dim=list(dimrivers, dimtime), missval=mv2, 
                      longname=paste("river total",cgem_vars[i],"load"))
  ncdata <- ncvar_add( ncdata, value )
}

for(i in 1:length(cgem_vars)){
# write the values
values <- array(0, dim=c(nr,nt))
for(itime in 1:nt){
  for(iriv in 1:nr){
    values[iriv,itime] <- cgem_init_values[i]
      #as.numeric(df[itime,iriv])
  }
}
ncvar_put(ncdata, varid=cgem_vars[i], vals=values, start=c(1,1), count=c(-1,-1), verbose=FALSE )
}

#Close netCDF file and clear data
nc_close(ncdata)
rm(ncdata)

#Old river file is 367.3KB and new one is 2.2MB