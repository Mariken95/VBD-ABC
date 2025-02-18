#________________________________________________________________________________
## Title:
# Functions for data preparation

## Description:
# This files contains functions that are used in data preparation steps 

## Functions
# 1. Create spatial dataframe
# 2. Add columns for mosquito abundance as placeholder
# 3. Estimate bird abundance
# 4. Extract temperature
# 5. Estimate mosquito abundance & emergence per time step


#_________________________________________________________________________________


# Blackbird abundance from Random Forest model ----------------------------

create.birdabundance <- function(abundancemaps, modelraster, outputfile, verbose) {
  if (verbose==1) { 
    browser()
  }
  
  modelraster <- rast(modelraster)
  
  for (i in 1:length(abundancemaps)) { 
    
  print(i)
    
  # 1. Load data
  abundancemap <- rast(abundancemaps[i])
  
  # 2. Turn raw bird unit raster into spatial dataframe to allow sum and mean to be calculated
  rasterBirdAb.DF <- as.data.frame(abundancemap, xy=TRUE) 

  # 3. Extract/transform bird units to fit raster used in model
  # crop bird units raster to be the same size as model raster
  rasterBirdAbCropped <- crop(abundancemap, modelraster)
  # resample to same grid as model raster
  rasterBirdAbresampled=resample(rasterBirdAbCropped,modelraster,method='bilinear') 

  
  # turn bird unit raster into spatial dataframe 
  croppedBirdAbraster.DF <- as.data.frame(rasterBirdAbresampled, xy=TRUE) 
  
  # calculate total number of bird units (sum of all raster values)
  print(sum(croppedBirdAbraster.DF[3], na.rm=TRUE)) 

  # 2879 bird units in total 
  # corresponds to 500,000-900,000 blackbirds https://www.vogelbescherming.nl/ontdek-vogels/kennis-over-vogels/vogelgids/vogel/merel
  # so a value of 1 unit in raster equals 700,000/2900 = 243 blackbirds 
  
  # 4. Turn bird unit values into real number of birds. see calculation above
  croppedBirdAbraster.DF[3] <- croppedBirdAbraster.DF[3]*243

  # 5. Add bird abundance to skeleton raster
  outputfile <- left_join(outputfile, croppedBirdAbraster.DF, by=c("x", "y")) 
  # if NA (bc bird raster misses some values on the border, slightly different shapefile),
  # take value from cell below
  birdabundance <- fill(outputfile, grep("ref_Merel", names(outputfile)), .direction="up")
  
  }
  
  return(birdabundance)
}

# Mosquito abundance from Random Forest model ----------------------------

create.mosquitoabundance <- function(abundancemaps, modelraster, outputfile, verbose=0) {
  if (verbose==1) { 
    browser()
  }
  
  for (i in 1:length(abundancemaps)) { #length(abundancemaps)
    
    print(i)
    
    # 1. Load data
    abundancemap <- raster(abundancemaps[i])

    # 2. Extract/transform to fit raster used in model
    # crop raster to be the same size as model raster
    rasterMosqAbCropped <- crop(abundancemap, extent(modelraster))
    # resample to same grid as model raster
    rasterMosqAbresampled=resample(rasterMosqAbCropped,modelraster,method='bilinear') 

    # 3. Multiply raw numbers for more realistic mosquito numbers
    # turn into spatial DF to allow joining
    croppedMosqAbraster.DF <- as.data.frame(rasterMosqAbresampled, xy=TRUE) 
    croppedMosqAbraster.DF[,-c(1:2)] <- croppedMosqAbraster.DF[,-c(1:2)] * 1000
    
    # 4. Add mosquito abundance to skeleton raster
    outputfile <- left_join(outputfile, croppedMosqAbraster.DF, by=c("x", "y")) 
    # if NA (bc mosquito raster misses some values on the border, slightly different shapefile),
    # take value from cell below
    mosquitoabundance <- fill(outputfile, grep("predictedr", names(outputfile)), .direction="up")
    
  }
  
  return(mosquitoabundance)
}


# Temperature -------------------------------------------------------------

# Loop through all files, return dataframe with temperature values (row=raster location, column=date)

create.temperature <- function(temperaturemaps, modelraster, outputfile) {
  
  for(i in 1:length(temperaturemaps)){ 
    print(i)
    
    tnc <- tidync::tidync(temperaturemaps[i])
    data <- hyper_tibble(tnc) %>%
      sf::st_as_sf(coords = c("x", "y"), crs = 28992) 
    
    # get x and y column instead of geometry list column (required to be able to use 'rasterFromXYZ' later)
    separated_coord <- data %>%
      mutate(x = unlist(map(data$geometry,1)),
             y = unlist(map(data$geometry,2)))
    
    # extract date from filename
    date <- as.character(substring(temperaturemaps[i],28,35))
    date <- as.Date(date, format="%Y%m%d")
    
    # only keep prediction and x, y columns
    data <- separated_coord[,c("x", "y", "prediction")] 
    data <- st_set_geometry(data, NULL) # also drop geometry column
    
    # extract/transform to fit raster
    # first turn both raster DF and temp DF into a raster object
    rasterTemp <- rasterFromXYZ(data)
    
    # crop temperature raster to be the same size as model raster
    rasterTempCropped <- crop(rasterTemp, extent(modelraster))
    # resample temperature to same grid as model raster
    rasterTempResampled=resample(rasterTempCropped,modelraster,method='bilinear')
    
    # add temperature values as column to spatialDF
    croppedTempraster.DF <- as.data.frame(rasterTempResampled, xy=TRUE) %>% drop_na(prediction)
    # rename prediction column to its date
    names(croppedTempraster.DF)[names(croppedTempraster.DF) == "prediction"] <- as.character(date)
    
    # add temperature to skeleton raster 
    outputfile <- left_join(outputfile, croppedTempraster.DF, by=c("x", "y")) 
  }  
  
  return(outputfile)
}


# Mosquito death & emergence ----------------------------------------------

# Death rate 

DeathRateCalc <- function(temperature) {
  ifelse(temperature <= 32, 1/(-4.86 * temperature + 169.8), 
         ifelse(temperature > 32, 0.17, 100)) # Shocket et al. cut-off own decision
}


create.mosquito_14d <- function(mosquitoabundance, timespan=92:274, verbose=0) {
  if (verbose==1) { 
    browser()
  }
  
  mosquitoabundance <- mosquitoabundance[,-c(1:3)]
  colnames(mosquitoabundance) <- paste(timespan)
  
    # take 14 day average to get more steady surface
    # first transpose so locations are in columns (take rolling mean over column)
  mosquitoabundanceTranspose <- as.data.frame(t(mosquitoabundance))
  mosquitoabundanceTranspose_14d <- rollapply(mosquitoabundanceTranspose, width=14, FUN=mean, fill="extend", by.column = TRUE)

    # transpose back so days are in columns (1 location is 1 row)
  mosquitoabundance_14d <- (t(mosquitoabundanceTranspose_14d))
  matplot(t(mosquitoabundance_14d[10:40,]),type="l")
  
  return(mosquitoabundance_14d)
}
 

create.mosquitoEmergence <- function(mosquitoabundance_14d, timespan=92:274, verbose=0) { 
  if (verbose==1) { 
    browser()
  }
  
  # Predefine outcome matrices 
  mosquitoesEmerging <- matrix(0,nrow(mosquitoabundance_14d),length(timespan)) # row = location, column = day. 

  # Derive emergence
  for (i in 1:nrow(mosquitoabundance_14d)){       # run model for each location
    DeltaN = diff(mosquitoabundance_14d[i,])
    
    for (j in 2:dim(mosquitoabundance_14d)[2]) {  # for each point in time, from day 2 because we can't derive delta for day 1
      Death = DeathRate[i,j]*mosquitoabundance_14d[i,j]
      Emergence  = DeltaN[j-1] + Death            # -1 due to the way diff is calculated (diff between col 1 and 2 is in position 1)
      mosquitoesEmerging[i,j] <- Emergence        # Number of mosquitoes to emerge 
      
    }
  }
  
  mosquitoesEmerging <- round(mosquitoesEmerging)
  matplot(t(mosquitoesEmerging[20:50,]), type="l")
  
  # Emergence can't be negative
  mosquitoesEmerging[mosquitoesEmerging<0] <- 0

  return(mosquitoesEmerging)
  
}


create.mosquitoDeath <- function(mosquitoabundance_14d, DeathRate) {
  mosquitoesDying <- mosquitoabundance_14d * DeathRate
  mosquitoesDying <- round(mosquitoesDying)

  return(mosquitoesDying)
  
}