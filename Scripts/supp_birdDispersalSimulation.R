## Simulate bird movement on a raster

library(tidyverse) 
library(sf)
library(raster)
library(sp)
library(stars)
library(spatialEco) #bearing.distance
library(tmap)   # static map of movement proportions

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Individual movement simulation ------------------------------------------

movementSimulation <- function(gridsizeM, nbirds, shape, scale, version, verbose=0) {
  if (verbose==1) browser()
  ## Create 5x5km grid of 100 squares in NL
  r <-
    raster(
      xmn = 125000,
      ymn = 425000,
      xmx = 130000,
      ymx = 430000,
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
  # mapview(ptsSpatial)
  
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
  # mapview(destSpatial)
  
  ## Determine % in each cell
  # First create large raster
  rlarge <-
    raster(
      xmn = 50000,
      ymn = 350000,
      xmx = 200000,
      ymx = 500000,
      resolution = gridsizeM)
  
  crs(rlarge) <- CRS("+init=EPSG:28992")
  rlarge <- setValues(rlarge, 1:ncell(rlarge)) # assign index number

  # Join raster with points.
  # Each point gets assigned which raster cell it occurs in
  rlargeStars <- stars::st_as_stars(rlarge)
  rlargePolygon <- st_as_sf(rlargeStars)
  st_crs(rlargePolygon) <- "+init=EPSG:28992" # assign crs
  
  pointsInRaster <- st_join(destSpatial, rlargePolygon)
  
  # Calculate proportion in each cell
  cellProp <- prop.table(table(pointsInRaster$layer))*100
  cellProp <- as.data.frame(cellProp) %>% rename(layer = Var1, proportion = Freq)
  
  # if visual output:
  cellProp$proportion <- round(cellProp$proportion, digits=1)

  # Merge with large raster to show % on map
  cellPropSpatial <- merge(rlargePolygon, cellProp, all.x = T)
  plotCellProp <- tm_shape(cellPropSpatial) + tm_polygons("proportion") +
    tm_text("proportion", size=0.6) + tm_layout(main.title = version, main.title.size=1)

 return(plotCellProp)
  
}



# Daily movement --------------------------------------------------------

# * baseline scenario -----------------------------------------------------

# within season - breeding: Weibull scale=329.4 & shape=0.44
visualmatrix.breeding <- movementSimulation(gridsizeM=5000, nbirds=1000000, shape=0.44, scale=329.4, version="Breeding") 

# within season - summer: Weibull scale=63.4 & shape=0.46
visualmatrix.summer <- movementSimulation(gridsizeM=5000, nbirds=1000000, shape=0.46, scale=63.4, version="Summer") 

tmap_arrange(visualmatrix.breeding, visualmatrix.summer, ncol=2, nrow= 1) 


# * scenario 20% ----------------------------------------------------------
# 20% of population is wild duck, 80% is blackbird

# breeding season
# wild duck scale parameter: 329.4*31.5 = 10376.1

# scale parameter 20% scenario:
0.8*329.4 + 0.2*10376.1 # 2338.74
median(rweibull(n=10000000, shape=0.44, scale=2338.74)) # 1016m
quantile(rweibull(n=10000000, shape=0.44, scale=2338.74), c(0.25, 0.5, 0.75)) 

# estimated median dispersal:
0.8*143+0.2*4500 # 1014m 
visualmatrix.br.sc20 <- movementSimulation(gridsizeM=5000, nbirds=100000, shape=0.44, scale=2338.74, version="Breeding sc 20%") 

# summer
# wild duck scale parameter: 63.4*52.5 = 3328.5

# scale parameter 20% scenario:
0.8*63.4 + 0.2*3328.5 # 716.42
median(rweibull(n=10000000, shape=0.46, scale=716.42)) # 323m
quantile(rweibull(n=10000000, shape=0.46, scale=716.42), c(0.25, 0.5, 0.75)) 
visualmatrix.su.sc20 <- movementSimulation(gridsizeM=5000, nbirds=100000, shape=0.46, scale=716.42, version="Summer sc 20%") 

# * plot comparison -------------------------------------------------------

# Breeding season
set.seed(123)  # Set seed for reproducibility
sample_bb <- data.frame(values = rweibull(1000, shape = 0.44, scale = 329.4))
sample_sc20 <- data.frame(values = rweibull(1000, shape = 0.44, scale = 2338.74))
breeding.season.samples <- cbind(sample_bb, sample_sc20, sample_bb) 

colnames(breeding.season.samples) <- c("Blackbird", "Increased \nmovement", "Sensitivity \nanalysis")
breeding.season.samples <- breeding.season.samples %>%
  pivot_longer(cols=c(1:3), names_to = "Scenario", values_to = "Value")
breeding.season.samples$Scenario <- factor(breeding.season.samples$Scenario, 
                                           levels = c("Blackbird", "Increased \nmovement", "Sensitivity \nanalysis"))

ggplot(data=breeding.season.samples, aes(x=Scenario, y=Value/1000, group=Scenario, fill=Scenario)) +
  geom_violin(scale="width") +
  ylim(0, 100) +
  labs(y = "Distance in KM", x="") +
  ggtitle("Breeding season") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=18), strip.text.x=element_text(size=18),
        plot.title=element_text(size=18), legend.title=element_blank(), legend.position = "none")
        
ggsave("../Output/Plots/dispersal/scenarios.breeding.png", 
       width=7, height=6)

# Summer season
set.seed(123)  # Set seed for reproducibility
sample_bb <- data.frame(values = rweibull(1000, shape = 0.46, scale = 63.4))
sample_sc20 <- data.frame(values = rweibull(1000, shape = 0.46, scale = 716.42))
sample_bb.sens <- data.frame(values = rweibull(1000, shape = 0.44, scale = 329.4))

summer.season.samples <- cbind(sample_bb, sample_sc20, sample_bb.sens) 
colnames(summer.season.samples) <- c("Blackbird", "Increased \nmovement", "Sensitivity \nanalysis")

summer.season.samples <- summer.season.samples %>%
  pivot_longer(cols=c(1:3), names_to = "Scenario", values_to = "Value")
summer.season.samples$Scenario <- factor(summer.season.samples$Scenario, 
                                         levels = c("Blackbird", "Increased \nmovement", "Sensitivity \nanalysis"))

ggplot(data=summer.season.samples, aes(x=Scenario, y=Value/1000, group=Scenario, fill=Scenario)) +
  geom_violin(scale="width") +
  ylim(0, 50) +
  labs(y = "Distance in KM", x="") +
  ggtitle("Post-breeding season") +
  scale_fill_brewer(palette = "Dark2") +
  theme_bw() +
  theme(axis.text=element_text(size=18), axis.title=element_text(size=18), strip.text.x=element_text(size=18),
        plot.title=element_text(size=18), legend.title=element_blank(), legend.position = "none")


ggsave("../Output/Plots/dispersal/scenarios.summer.png", 
       width=7, height=6)

# Seasonal movement -------------------------------------------------------

# natal dispersal - Weibull scale=668.2, shape=0.64
visualmatrix.seasonal.juv <- movementSimulation(gridsizeM=5000, nbirds=1000000, shape=0.64, scale=668.2, version="Natal dispersal") # 86% remains, whole grid filled

# breeding dispersal - Weibull scale=459.4, shape=0.39
visualmatrix.seasonal.adu <- movementSimulation(gridsizeM=5000, nbirds=1000000, shape=0.39, scale=459.4, version="Breeding dispersal") # 93% remains, whole grid filled

tmap_arrange(visualmatrix.seasonal.juv, visualmatrix.seasonal.adu, ncol=2, nrow= 1) 
