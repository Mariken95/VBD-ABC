#________________________________________________________________________________
## Title:
# createRaster

## Description:
# This files creates a spatial raster for the analysis/paper called: 

#________________________________________________________________________________

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tmap)       # provides NLD_muni and NLD_prov maps (shapefiles). does not contain water-areas
library(tidyverse) 
library(Rcpp)
library(sf)
library(raster)
library(rgdal)
library(viridis)    # fill
library(cowplot)    # plot grids
library(mapview)

# Create spatial DF without water regions --------------------------------------------------------
# in this version we assigned index numbers to the whole square. so we could extract which numbers are 
# outside of the square. these cols and rows can then be deleted from the ldata dispersal df

# define size of raster
km=5

create.raster.nowater <- function(gridsize.km, verbose=0) {
  if (verbose==1) { 
    browser()
  }
  
  # create empty raster
  r=raster()
  r <- raster(nrows=68, ncols=56, xmn= 646, ymn= 296050, xmx = 278975, ymx = 636456, resolution = gridsize.km*1000)
  
  # set crs
  crs(r)<- CRS("+init=EPSG:28992")
  r <- setValues(r, seq(from=1, to=ncell(r), by=1)) 
  
  # get NL shape and only keep these cells 
  data("NLD_prov")
  
  # turn country map into raster. values outside country are NA
  nl_r = rasterize(NLD_prov, r, getCover=TRUE)
  nl_r[nl_r<0.5] <- NA # only keep cells that are for more than x% (used 50% for 5x5km raster) within the country
  
  # for original raster, hide the squares that are outside the country 
  r_masked=mask(r, nl_r)
  
  # crop r_masked according to NL boundaries
  croppedNLraster <- crop(x=r_masked, y=extent(NLD_prov))

  # save raster as spatial data frame
  croppedNLraster.DF <- as.data.frame(croppedNLraster, xy=TRUE) %>% drop_na(layer)
  # assign index to each grid cell 
  croppedNLraster.DF <- croppedNLraster.DF %>% mutate(index=seq(from=1, to=nrow(croppedNLraster.DF)))
  
  # pick relevant output format
  return(croppedNLraster)
  return(croppedNLraster.DF)
} 

croppedNLraster <- create.raster.nowater(gridsize.km = 5)

write.csv(croppedNLraster.DF,"C:/Repos/vbd-siminf/Data/spatialDF5kmNoWater.csv", row.names = FALSE)
writeRaster(croppedNLraster, "C:/Repos/vbd-siminf/Data/spatialDF5kmNoWater", format="GTiff", overwrite=TRUE)

