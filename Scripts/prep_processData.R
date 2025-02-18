#________________________________________________________________________________
## Title:
# Data Processing 

## Description:
# This files processes intput data so it can be fed into the model 

# 1. Create files for starting conditions
# 2. Create events data frames
# 3. Create matrix of dispersal probabilities

## Author:
# Mariken de Wit
# mariken.dewit@wur.nl

## Date: 1-12-2022

#_________________________________________________________________________________

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(sf)
library(raster)
library(sp)
library(stars)
library(spatialEco) # for bearing.distance

source('functions_processData.R')


# Create starting abundance files -----------------------------------------

load("../Data/model_input/birdabundance5kmNoWater.RData")
load("../Data/model_input/mosquitoabundance_14d_5kmNoWater_2016.RData")

startingabundance <- create.startingabundance(birdabundance = birdabundance, 
                                              mosquitoabundance = mosquitoabundance_14d, 
                                              nrDatasets=1, verbose=0)

names(startingabundance)[1] <- "mean" # first list is mean
save(startingabundance, file = "../Data/model_input/startingabundance.RData")




# DISPERSAL ---------------------------------------------------------------

# see data_birdDispersalSimulation.R for viz and choice of dispersal parameters

# Use distance between raster cells

# Approach
# 1. create matrix with values of distance between raster cells
# 2. replace distance values with probability to move that distance

modelrasterDF <- read.csv("../Data/model_input/spatialDF5kmNoWater.csv")
modelrasterDF$node <- 1:nrow(modelrasterDF)

# creates dataframe from all combinations of supplied vectors
distDF=expand.grid(source=modelrasterDF$node,dest=modelrasterDF$node) 

# add coordinates source
distDF = modelrasterDF %>% 
  mutate(source=node, source_x=x, source_y=y) %>%
  dplyr::select(source, source_x, source_y) %>% left_join(distDF,.,by="source")

# add coordinates dist
distDF = modelrasterDF %>% 
  mutate(dest=node, dest_x=x, dest_y=y) %>%
  dplyr::select(dest, dest_x,dest_y) %>% left_join(distDF,.,by="dest")

# calculate distances
distDF=distDF %>%
  rowwise() %>% 
  mutate(distance=getDistance(c(dest_x,dest_y),c(source_x,source_y), resolution=1))

# Dispersal: DAILY --------------------------------------------------------

# * Scenario: Baseline ----------------------------------------------------

# ** breeding season ----

mov.breeding.baseline <- movementSimulation_withDist(gridsizeM=5000, nbirds=1000000, shape=0.44, scale=329.4, verbose=0) 
matrix.breeding.baseline <- create.ldatamatrix(distDF = distDF, movement.output=mov.breeding.baseline, n.nodes=1398)

# set column names
colnames.breeding <- c(1:ncol(matrix.breeding.baseline))
colnames.breeding <- sub("^","br",colnames.breeding)
colnames(matrix.breeding.baseline) <- colnames.breeding

# ** summer season ----

mov.summer.baseline <- movementSimulation_withDist(gridsizeM=5000, nbirds=1000000, shape=0.46, scale=63.4, verbose=0) 
matrix.summer.baseline <- create.ldatamatrix(distDF = distDF, movement.output=mov.summer.baseline, n.nodes=1398)

# set column names
colnames.summer <- c(1:ncol(matrix.summer.baseline))
colnames.summer<-sub("^","sum",colnames.summer)
colnames(matrix.summer.baseline) <- colnames.summer

# save each scenario separately
# blackbird scenario
dispersal_daily.bb <- list(matrix.breeding.baseline, matrix.summer.baseline)
names(dispersal_daily.bb)[[1]] <- "breeding"
names(dispersal_daily.bb)[[2]] <- "summer"

save(dispersal_daily.bb, file="../Data/model_input/dispersal_daily.bb.RData")


# * Scenario: 20% other birds ---------------------------------------------------

# ** breeding season ----

mov.breeding.sc20 <- movementSimulation_withDist(gridsizeM=5000, nbirds=1000000, shape=0.44, scale=2338.74, verbose=0) 
matrix.breeding.sc20 <- create.ldatamatrix(distDF = distDF, movement.output=mov.breeding.sc20, n.nodes=1398)

# set column names
colnames(matrix.breeding.sc20) <- colnames.breeding

# ** summer season ----

mov.summer.sc20 <- movementSimulation_withDist(gridsizeM=5000, nbirds=1000000, shape=0.46, scale=716.42, verbose=0) 
matrix.summer.sc20 <- create.ldatamatrix(distDF = distDF, movement.output=mov.summer.sc20, n.nodes=1398)

# set column names
colnames(matrix.summer.sc20) <- colnames.summer

# save each scenario separately
# 20% other bird scenario
dispersal_daily.sc20 <- list(matrix.breeding.sc20, matrix.summer.sc20)
names(dispersal_daily.sc20)[[1]] <- "breeding"
names(dispersal_daily.sc20)[[2]] <- "summer"

save(dispersal_daily.sc20, file="../Data/model_input/dispersal_daily.sc20.RData")


# * Update colnames -------------------------------------------------------
# update column names to contain scenario name to avoid duplicate column names in ldata 
# required when using dispersal scenarios for multi-host system
load("../Data/model_input/dispersal_daily.bb.RData")
old_colnames <- colnames(dispersal_daily.bb[["breeding"]])
new_colnames <- paste0("bb_", old_colnames)
colnames(dispersal_daily.bb[["breeding"]]) <- new_colnames

old_colnames <- colnames(dispersal_daily.bb[["summer"]])
new_colnames <- paste0("bb_", old_colnames)
colnames(dispersal_daily.bb[["summer"]]) <- new_colnames

save(dispersal_daily.bb, file="../Data/model_input/dispersal_daily.bb.RData")


load("../Data/model_input/dispersal_daily.sc20.RData")
old_colnames <- colnames(dispersal_daily.sc20[["breeding"]])
new_colnames <- paste0("sc20_", old_colnames)
colnames(dispersal_daily.sc20[["breeding"]]) <- new_colnames

old_colnames <- colnames(dispersal_daily.sc20[["summer"]])
new_colnames <- paste0("sc20_", old_colnames)
colnames(dispersal_daily.sc20[["summer"]]) <- new_colnames

save(dispersal_daily.sc20, file="../Data/model_input/dispersal_daily.sc20.RData")


# Dispersal: BETWEEN SEASON --------------------------------------------------------

# * Scenario: Baseline -------------------------------------------------------------

# ** juvenile ----

mov.juv.baseline <- movementSimulation_withDist(gridsizeM=5000, nbirds=1000000, shape=0.64, scale=668.2, verbose=0) 
seasonal.juv.baseline <- create.seasonalDF(distDF = distDF, movement.output=mov.juv.baseline)

# ** adult --------

mov.adu.baseline <- movementSimulation_withDist(gridsizeM=5000, nbirds=1000000, shape=0.39, scale=459.4, verbose=0) 
seasonal.adu.baseline <- create.seasonalDF(distDF = distDF, movement.output=mov.adu.baseline)

dispersal_seasonal.baseline <- list(seasonal.juv.baseline, seasonal.adu.baseline)
names(dispersal_seasonal.baseline)[[1]] <- "juvenile"
names(dispersal_seasonal.baseline)[[2]] <- "adult"
save(dispersal_seasonal.baseline, file="../Data/model_input/dispersal_seasonal.bb.RData")

# Events ------------------------------------------------------------------

# Each dispersal scenario has its own events dataframe. 
# These are saved separately with the scenario indicated in the name

# * Mosquito emergence -------------------------

eventsMosqEmerg <- list()
temp <- list()

year <- c(2016, 2017, 2018, 2019, 2020, 2021, 2022)

load("../Data/model_input/mosquitoesEmerging_5kmNoWater_2016.RData")

for (i in 1:length(mosquitoesEmerging)) {
  
  for (year in year) {
    print(year)
    
    load(paste0("../Data/model_input/mosquitoesEmerging_5kmNoWater_", year, ".RData"))
    
    emergenceDF <- as.data.frame(mosquitoesEmerging[[i]])
    nodes <- 1:nrow(emergenceDF)
    
    timespan <- 1:ncol(emergenceDF)
    timespan <- timespan + 91 # so that the first day is day 92 = 1 april 2016
    
    if (year == 2017) {
      timespan <- timespan + 365 
    } else if (year == 2018) {
      timespan <- timespan + 365*2
    } else if (year == 2019) {
      timespan <- timespan + 365*3
    } else if (year == 2020) {
      timespan <- timespan + 365*4+1
    } else if (year == 2021) {
      timespan <- timespan + 365*5+1
    } else if (year == 2022) {
      timespan <- timespan + 365*6+1
    }
    
    print(timespan)
    
    firstday <- timespan[1]-1 #add new mosquitoes at the last day before new season starts
    
    emergence.unlist <- unlist(emergenceDF)
    
    eventsMosqEmerg.i <- data.frame(event      = rep("enter", length(nodes)*length(timespan)),
                                    time       = rep(timespan, each = length(nodes)),
                                    node       = rep(nodes, length(timespan)),
                                    dest       = 0,
                                    n          = emergence.unlist, 
                                    proportion = 0,
                                    select     = 1,
                                    shift      = 0)
    
    if (year == 2016) {
      temp <- rbind(temp, eventsMosqEmerg.i)
    }
    
    if (year != 2016) {
      
      # on first day of new year, number of mosquito needs to equal model predictions 
      load(paste0("../Data/model_input/mosquitoabundance5kmNoWater_", year, ".RData"))
      
      eventsMosqEmerg.firstday <- data.frame(event      = rep("enter", length(nodes)),
                                             time       = rep(firstday, length(nodes)),
                                             node       = nodes,
                                             dest       = 0,
                                             n          = round(unlist(mosquitoabundance[[i]][4])), 
                                             proportion = 0,
                                             select     = 1,
                                             shift      = 0)
      total.eventsMosqEmerg <- list()
      total.eventsMosqEmerg <- rbind(eventsMosqEmerg.i, eventsMosqEmerg.firstday)
      
      temp <- rbind(temp, total.eventsMosqEmerg)
      
    }
    
    eventsMosqEmerg[[i]] <- temp
    
  }
}

names(eventsMosqEmerg)[1] <- "mean" # first list is mean


# * Mosquito death -----------------------------

eventsMosqDeath <- list()
temp <- list()

year <- c(2016, 2017, 2018, 2019, 2020, 2021, 2022)

load("../Data/model_input/mosquitoesDying_5kmNoWater_2016.RData")

for (year in year) {
  print(year)
  
  load(paste0("../Data/model_input/mosquitoesDying_5kmNoWater_", year, ".RData"))
  
  for (i in 1:length(mosquitoesDying)) {
    
    deathDF <- as.data.frame(mosquitoesDying[[i]])
    nodes <- 1:nrow(deathDF)
    timespan <- 1:ncol(deathDF)
    timespan <- timespan + 91 # so that the first day is day 92 = 1 april 2016
    diapauseEnd <- 289 # 15th of October
    
    
    if (year == 2017) {
      timespan <- timespan + 365 
      diapauseEnd <- diapauseEnd + 365
    } else if (year == 2018) {
      timespan <- timespan + 365*2
      diapauseEnd <- diapauseEnd + 365*2
    } else if (year == 2019) {
      timespan <- timespan + 365*3
      diapauseEnd <- diapauseEnd + 365*3
    } else if (year == 2020) {
      timespan <- timespan + 365*4+1
      diapauseEnd <- diapauseEnd + 365*4+1
    } else if (year == 2021) {
      timespan <- timespan + 365*5+1
      diapauseEnd <- diapauseEnd + 365*5+1
    } else if (year == 2022) {
      timespan <- timespan + 365*6+1
      diapauseEnd <- diapauseEnd + 365*6+1
    }  
    
    death.unlist <- unlist(deathDF)
    
    eventsMosqDeath.i <- data.frame(event      = rep("exit", length(nodes)*length(timespan)),
                                    time       = rep(timespan, each = length(nodes)),
                                    node       = rep(nodes, length(timespan)),
                                    dest       = 0,
                                    n          = death.unlist,
                                    proportion = 0,
                                    select     = 2,
                                    shift      = 0)
    
    eventsMosqDeath.diapauseS <- data.frame(event      = rep("exit", length(nodes)),
                                            time       = rep(diapauseEnd, length(nodes)),
                                            node       = nodes,
                                            dest       = 0,
                                            n          = 0,
                                            proportion = 0.99,
                                            select     = 1,
                                            shift      = 0)
    
    eventsMosqDeath.diapauseE <- data.frame(event      = rep("exit", length(nodes)),
                                            time       = rep(diapauseEnd, length(nodes)),
                                            node       = nodes,
                                            dest       = 0,
                                            n          = 0,
                                            proportion = 1,
                                            select     = 3,
                                            shift      = 0)
    
    eventsMosqDeath.diapauseI <- data.frame(event      = rep("exit", length(nodes)),
                                            time       = rep(diapauseEnd, length(nodes)),
                                            node       = nodes,
                                            dest       = 0,
                                            n          = 0,
                                            proportion = 1,
                                            select     = 4,
                                            shift      = 0)
    
    
    total.eventsMosqDeath <- list()
    total.eventsMosqDeath <- rbind(eventsMosqDeath.i, eventsMosqDeath.diapauseS, eventsMosqDeath.diapauseE, eventsMosqDeath.diapauseI)
    
    temp <- rbind(temp, total.eventsMosqDeath)
    
    eventsMosqDeath[[i]] <- temp
    
    
  }
}
names(eventsMosqDeath)[1] <- "mean" # first list is mean


# * Ageing ----------------------------------------------------------------
# on 30th of April, all juvenile birds become adults (births start 1st of May)

# day number of 30th of April
endApril <- c(486, 486+365, 486+365*2, 486+365*3+1, 486+365*4+1, 486+365*5+1)
nodes <- 1:nrow(emergenceDF)

eventsAgeing.bb <- data.frame(event = rep("intTrans", length(nodes)*length(endApril)),
                              time  = rep(endApril, each = length(nodes)),
                              node  = rep(nodes, length(endApril)),
                              dest  = 0,
                              n     = 0, 
                              proportion = 1,
                              select = 6,
                              shift  = 1)

eventsAgeing.res <- data.frame(event = rep("intTrans", length(nodes)*length(endApril)),
                               time  = rep(endApril, each = length(nodes)),
                               node  = rep(nodes, length(endApril)),
                               dest  = 0,
                               n     = 0, 
                               proportion = 1,
                               select = 10,
                               shift  = 4)

# * Seasonal dispersal ----------------------------------------------------------------
# on 1st of July, seasonal/natal dispersal takes place for all birds 

load("../Data/model_input/dispersal_seasonal.bb.RData")
dispersal_seasonal <- dispersal_seasonal.baseline


# day number of 1st of July
dispersalDay <- c(183, 183+365, 183+365*2, 183+365*3+1, 183+365*4+1, 183+365*5+1, 183+365*6+1)

eventsDispersal_juv <- data.frame(event = "extTrans",
                                  time  = rep(dispersalDay, each = nrow(dispersal_seasonal[[1]])),
                                  node  = rep(dispersal_seasonal[[1]]$source, length(dispersalDay)),
                                  dest  = rep(dispersal_seasonal[[1]]$dest, length(dispersalDay)),
                                  n     = 0, 
                                  proportion = rep(dispersal_seasonal[[1]]$proportion, length(dispersalDay)),
                                  select = 6,
                                  shift  = 2)

eventsDispersal_adu <- data.frame(event = "extTrans",
                                  time  = rep(dispersalDay, each = nrow(dispersal_seasonal[[2]])),
                                  node  = rep(dispersal_seasonal[[2]]$source, length(dispersalDay)),
                                  dest  = rep(dispersal_seasonal[[2]]$dest, length(dispersalDay)),
                                  n     = 0, 
                                  proportion = rep(dispersal_seasonal[[2]]$proportion, length(dispersalDay)),
                                  select = 7,
                                  shift  = 2)

eventsDispersal <- rbind(eventsDispersal_juv, eventsDispersal_adu)

# * Events data frame ------------------------

events <- list()
for (i in 1:length(eventsMosqEmerg)) {
  events[[i]] <- rbind(eventsMosqEmerg[[i]], eventsMosqDeath[[i]], eventsAgeing.bb, eventsAgeing.res, eventsDispersal)
  
}

events_dispersal.bb <- events
save(events_dispersal.bb, file="../Data/model_input/events_5kmNoWater_2016_22_agestruc_dispersal.bb.RData")
