#________________________________________________________________________________
## Title:
# Functions for preparing data to feed into model 

## Description:
# Functions to create abundance starting conditions (u0), events, and dispersal probabilities

#_________________________________________________________________________________



# Abundance starting conditions -------------------------------------------

# starting conditions should be a list with each dataframe starting with x, y, layer

create.startingabundance <- function(birdabundance, mosquitoabundance, nrDatasets, verbose=0) {
  
  if (verbose==1) { 
    browser()
  }
  
  # read in bird abundance file
  # select first data column
  # read in mosquito abundance files
  # select first data column from each file
  
  startingabundance <- list() # Create list to store all starting abundance dataframes
  
  # loop by selecting next column in bird abundance file and next file in mosquito abundance files
  
  for (i in 1:nrDatasets) {
    startingabundance[[i]] <- birdabundance[1:3]    # add x, y, layer column
    startingabundance[[i]][4] <- birdabundance[i+3] # add abundance data to 4th column
    colnames(startingabundance[[i]])[4] <- "blackbird_abundance"
    
    mosquitoabundance.i <- mosquitoabundance[[i]] # select file
    mosquitoabundance.i <- mosquitoabundance.i[,1]       # select column with abundance estimates at t=0
    startingabundance[[i]][5] <- mosquitoabundance.i
    colnames(startingabundance[[i]])[5] <- "mosquito_abundance"
    
  }
  
  return(startingabundance)
  
}



# Calculate distance ---------------------------------------------------------------

# function to calculate the distance between two cells, based on their coordinates
# used for dispersal probabilities

getDistance <- function (point.a, point.b, resolution = 1){
  distance.ab <- 0
  dx <- point.a[1] - point.b[1]
  dy <- point.a[2] - point.b[2]
  xdist <- dx * resolution
  ydist <- dy * resolution
  distance.ab <- sqrt(xdist^2 + ydist^2)
  return(distance.ab)
}


# Calculate % moved to each cell ------------------------------------------

# function that calculates what % of birds end up where.
# this function works on model raster so that it can be combined with calculated distances between nodes
# so that the % can immediately be merged based on distance

movementSimulation_withDist <- function(gridsizeM, nbirds, shape, scale, version, verbose=0) {
  ## Create 1 5x5km grid cell in NL 
  r <-
    raster(
      xmn = 145646,
      ymn = 471456-5000,
      xmx = 145646+5000,
      ymx = 471456,
      resolution = gridsizeM
    )
  crs(r) <- CRS("+init=EPSG:28992")
  r <- setValues(r, 1:ncell(r)) # assign index number
  
  ## Create points within raster
  p <- as(r@extent, 'SpatialPolygons')
  pts <- spsample(p, n = nbirds, "random")
  pts <- as.data.frame(pts)
  
  # convert to spatial DF
  ptsSpatial <- st_as_sf(pts, coords = c("x", "y"))
  st_crs(ptsSpatial) <- "+init=EPSG:28992" # assign crs

  ## Move individuals
  # draw from parameterised weibull distribution
  dest <-
    bearing.distance(pts$x, pts$y,
                     rweibull(nrow(pts), shape = shape, scale = scale),
                     runif(nrow(pts), min = 0, max = 360))
  
  dest <- as.data.frame(dest)
  dest <- dest %>% mutate(birdID = c(1:nrow(dest)))
  
  destSpatial <- st_as_sf(dest, coords=c("X", "Y"))
  st_crs(destSpatial) <- "+init=EPSG:28992" # assign crs

  ## Determine % in each cell
  # First load model raster of NL
  rlarge <- raster("../Data/model_input/spatialDF5kmNoWater.tif")
  
  crs(rlarge) <- CRS("+init=EPSG:28992")
  rlarge[!is.na(rlarge[])] <- c(1:length(rlarge[!is.na(rlarge[])])) # only set values for grid cells within NL surface (non-NA)

  # Join NL raster with points of where birds were moved to.
  # Each point gets assigned which raster cell it occurs in
  rlargeStars <- stars::st_as_stars(rlarge)
  rlargePolygon <- st_as_sf(rlargeStars)
  st_crs(rlargePolygon) <- "+init=EPSG:28992" # assign crs
  
  if (verbose==1) browser()
  
  pointsInRaster <- st_join(destSpatial, rlargePolygon)
  
  # Calculate proportion in each cell
  cellProp <- pointsInRaster %>% group_by(layer) %>% summarise(total = n()) %>% 
    mutate(proportion = total/nbirds) %>%
    dplyr::select(layer, proportion)
  
  cellProp <- st_set_geometry(cellProp, NULL)
  
  return(cellProp)
  
}



# Use % moved for daily dispersal: ldata ----------------------------------


# function to go from simulated dispersal to ldata format of movement matrix 

create.ldatamatrix <- function(distDF, movement.output, n.nodes) {
  
  # 1. combine the % that moves to each cell with calculated distances between cells

  dist_prop <- distDF %>% filter(source==694) # 694 is index nr of node where all birds started before they moved
  
        # add proportion to distDF (creates link between proportion & distance based on destination cell)
  dist_prop <- left_join(dist_prop, movement.output, by=c("dest"= "layer"))
  
        # keep unique distance values so we get 1 proportion per distance (avoid small inconsistencies from simulating)
        # remove % that stays, this will be added later to make sum equal to 1
  dist_prop_distinct <- dist_prop %>% filter(!is.na(proportion)) %>% 
    distinct(distance, .keep_all = T) %>%
    dplyr::select(distance, proportion) %>%
    filter(distance!=0)       
  
        # merge df with unique combinations into full df
  distDF.withProp <- left_join(distDF, dist_prop_distinct, by="distance")
  

  
  # 2. organise into matrix for ldata (all locations on x and y, cell = dispersal %)
  matrix.prop <- distDF.withProp[1:n.nodes,]
  
  for (i in 2:nrow(matrix.prop)) {
    
    newcol <- distDF.withProp$proportion[((i-1)*n.nodes+1):(n.nodes*i)]
    matrix.prop <- cbind(matrix.prop, newcol)
  }
  
  matrix.prop <- matrix.prop[,-c(1:7)]
  matrix.prop <- matrix.prop %>% replace(is.na(.),0)
  
  
  
  # 3. make sure that all columns add up to 1, % that stays is 1-prop.leaving
  for (i in 1:ncol(matrix.prop)) {
    
    matrix.prop[i,i] = 1-sum(matrix.prop[,i])
  }
  
  return(matrix.prop)
}

# Use % moved for seasonal dispersal: events ----------------------------------

create.seasonalDF <- function(distDF, movement.output) {
  
  # 1. combine the % that moves to each cell with calculated distances between cells
  
  dist_prop <- distDF %>% filter(source==694) # 694 is index nr of node where all birds started before they moved
  
  # add proportion to distDF (creates link between proportion & distance based on destination cell)
  dist_prop <- left_join(dist_prop, movement.output, by=c("dest"= "layer"))
  
  # keep unique distance values so we get 1 proportion per distance (avoid small inconsistencies from simulating)
  # remove % that stays, this will be added later to make sum equal to 1
  dist_prop_distinct <- dist_prop %>% filter(!is.na(proportion)) %>% 
    distinct(distance, .keep_all = T) %>%
    dplyr::select(distance, proportion) %>%
    filter(distance!=0)       
  
  # merge df with unique combinations into full df
  distDF.withProp <- left_join(distDF, dist_prop_distinct, by="distance")
  
  # remove rows where source == destination & with proportion < 0.1% (=0.001) & with NA proportion
  distDF.withProp_filtered <- distDF.withProp %>% 
    filter(source!=dest) %>%
    filter(proportion>0.001) %>%
    filter(!is.na(proportion))
  
  return(distDF.withProp_filtered)
}
