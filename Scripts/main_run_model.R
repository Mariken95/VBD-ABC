#________________________________________________________________________________

## Title:
# functions_main_Simulation

## Description:
# This files runs the simulation model for the analysis/paper called: 

#________________________________________________________________________________


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(SimInf)
library(tidyverse)

# set-up ------------------------------------------------------------------
rm(list=ls())

load("../Data/model_input/events_5kmNoWater_2016_22_agestruc_dispersal.bb.RData") # events
load("../Data/model_input/startingabundance.RData") # starting abundance


# SELECT RELEVANT FILES TO SOURCE DEPENDING ON MODEL VERSION

# * one-host model --------------------------------------------------------
source("prep_parameters_onehost.R")      # loads compartments, transitions, parameters, initial prevalence, ldata and v0
source("functions_run_model_onehost.R")   # loads 'buildModel function'


# two-host model ----------------------------------------------------------
source("prep_parameters_twohost.R")  # loads compartments, transitions, parameters, initial prevalence, ldata and v0
source("functions_run_model_twohost.R")           # loads 'buildModel function'


# set parameters ----------------------------------------------------------

# specify which dispersal scenario and therefore which ldata and events version to use
scenario.dispersal <- "blackbird"
# scenario.dispersal <- "sc20"
# scenario.dispersal <- "sensitivity"

abundanceReservoir <- "blackbird"
introduction <- "south"


if (scenario.dispersal == "blackbird") {
  ldata <- ldata_dispersal.bb
  events <- events_dispersal.bb
  rm(ldata_dispersal.bb, ldata_dispersal.sc20, ldata_dispersal.sens, events_dispersal.bb)
  
} else if (scenario.dispersal == "sc20") {
  ldata <- ldata_dispersal.sc20
  events <- events_dispersal.bb
  rm(ldata_dispersal.bb, ldata_dispersal.sc20, ldata_dispersal.sens, events_dispersal.bb)
  
} else if (scenario.dispersal == "sensitivity") {
  ldata <- ldata_dispersal.sens
  events <- events_dispersal.bb
  rm(ldata_dispersal.bb, ldata_dispersal.sc20, ldata_dispersal.sens, events_dispersal.bb)
  
} else {
  print("no valid dispersal scenario specified")
}

# set parameters for identifiability datasets
ldata$vertical_transm <- 0.1
ldata$transmissionProbHM <- 0.6
scalingParameter <- 0.1
ldata$disIndMortalityH <- 0.4
FOI.n <- rep(0.01, n/3)
FOI.m <- rep(0.01, n/3)
FOI.s <- rep(0.1, n/3)
FOI <- c(FOI.n, FOI.m, FOI.s)

# Run one-host model ---------------------------------------------------------------

## remove the events related to reservoir ageing (using select value 10) to be able to use the events in the one-host system
events[[1]] <- events[[1]][events[[1]]$select < 8, ]

model <- buildModel(transitions = transitions, compartments = compartments, ldata = ldata, E = E, N = N, v0 = v0, tspan = 92:2550, 
                    startingabundance = startingabundance, introduction = introduction, scalingParameter = scalingParameter, events = events, 
                    pts_fun = pts_fun, age_structure="yes", runtype ="mean", FOI=FOI, useFOI = "yes", verbose=0) 

# Run two-host model ---------------------------------------------------------------

model <- buildModel(transitions = transitions, compartments = compartments, ldata = ldata, E = E, N = N, v0 = v0, tspan = 92:2550, 
                    startingabundance = startingabundance, introduction = introduction, scalingParameter = scalingParameter, events = events, 
                    pts_fun = pts_fun, age_structure="yes", runtype ="mean", FOI=FOI, useFOI = "yes", 
                    abundanceReservoir = abundanceReservoir, verbose=0) 

results <- run(model)
traject <- trajectory(results)

save(traject, file="../Data/abc_identifiability/traject_E.Rda")



