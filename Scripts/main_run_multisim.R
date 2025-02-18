
rm(list=ls())

# test 
print('is something happening?')

.libPaths( c("/home/mdwit/R/library" , .libPaths() ) )

library(tidyverse,lib="/home/mdwit/R/library")
library(SimInf,lib="/home/mdwit/R/library")

# set-up ------------------------------------------------------------------

load("../Data/model_input/events_5kmNoWater_2016_22_agestruc_dispersal.bb.RData") # events
load("../Data/model_input/startingabundance.RData") # starting abundance


# set parameters ----------------------------------------------------------

args=(commandArgs(TRUE))
print(args)
if (length (args) == 0) {
  print("No arguments supplied.")
} else {
  for (i in 1:length(args)) {
    eval (parse (text = args[[i]] ))
  }
}

print(nr.runs)
print(scenario.dispersal)
print(scenario.parameters)
print(posterior)
print(outputPath)


# load posterior
load(paste(outputPath, "/", posterior, ".Rda", sep=""))

# select one-host or two-host model

if (scenario.parameters == "A") {
  # * one-host model
  source("prep_parameters_onehost.R")       # loads compartments, transitions, parameters, initial prevalence, ldata and v0
  source("functions_run_model_onehost.R")   # loads 'buildModel function'

  } else {
  # two-host model
  source("prep_parameters_twohost.R")         # loads compartments, transitions, parameters, initial prevalence, ldata and v0
  source("functions_run_model_twohost.R")     # loads 'buildModel function'
  
}

n <- nrow(ldata_dispersal.bb) # no of locations


# Create output filenames --------------------------------------------------------

# simulation output
filename = sprintf(paste(outputPath, 'multirun_%s.Rdata', sep=''),indexnr)

# Select scenario input ---------------------------------------------------------------

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
  stop()
}

# Sample from posterior & run model ---------------------------------------------------

selection <- round(runif(nr.runs, min=0, max=nrow(posterior.allchains)))
posteriorSamples <- as.data.frame(posterior.allchains[selection,])

listofresults <- list() # Create list to store all results dataframes

for(ii in 1:nr.runs) {
  
  print(ii)
  

# * specify parameters ----------------------------------------------------

  if(scenario.parameters == "A") {
    
    introduction <- "south"
    
    ldata$vertical_transm <- rep(posteriorSamples[ii, 'vertical_transm'], n)
    ldata$transmissionProbHM <- rep(posteriorSamples[ii, 'transmissionProbHM'], n)
    scalingParameter <- posteriorSamples[ii,'scalingParameter']
    ldata$disIndMortalityH <- rep(posteriorSamples[ii,'disIndMortality'], n)
    FOI <- rep(posteriorSamples[ii, 'FOI'], n)
    
  } else if(scenario.parameters == "B") {
    
    abundanceReservoir <- "blackbird"
    introduction <- "south"
    
    ldata$vertical_transm <- rep(posteriorSamples[ii, 'vertical_transm'], n)
    ldata$transmissionProbHM <- rep(posteriorSamples[ii, 'transmissionProbHM'], n)
    scalingParameter <- posteriorSamples[ii,'scalingParameter']
    ldata$disIndMortalityH <- rep(posteriorSamples[ii,'disIndMortality'], n)
    FOI <- rep(posteriorSamples[ii, 'FOI'], n)
    ldata$pref <- rep(posteriorSamples[ii,'bitingPref'], n)
    ldata$disIndMortalityR <- rep(posteriorSamples[ii,'disIndMortality'], n)  # reservoir disease-induced mortality equal to blackbird
    ldata$lifespanR <- 1/ldata$mortalityRateHAdu                              # reservoir adult mortality rate equal to blackbird
    
  } else if(scenario.parameters == "C") {
    
    abundanceReservoir <- "blackbird"
    introduction <- "south"
    
    ldata$vertical_transm <- rep(posteriorSamples[ii, 'vertical_transm'], n)
    ldata$transmissionProbHM <- rep(posteriorSamples[ii, 'transmissionProbHM'], n)
    scalingParameter <- posteriorSamples[ii,'scalingParameter']
    ldata$disIndMortalityH <- rep(posteriorSamples[ii,'disIndMortality'], n)
    FOI <- rep(posteriorSamples[ii, 'FOI'], n)
    ldata$lifespanR <- rep(posteriorSamples[ii,'lifespanR'], n)
    ldata$pref <- rep(posteriorSamples[ii,'bitingPref'], n)
    ldata$disIndMortalityR <- 0
    
  } else if(scenario.parameters == "D") {
    
    abundanceReservoir <- "blackbird"
    introduction <- "south"
    
    ldata$vertical_transm <- rep(posteriorSamples[ii, 'vertical_transm'], n)
    ldata$transmissionProbHM <- rep(posteriorSamples[ii, 'transmissionProbHM'], n)
    scalingParameter <- posteriorSamples[ii,'scalingParameter']
    ldata$disIndMortalityH <- rep(posteriorSamples[ii,'disIndMortality'], n)
    FOI <- rep(posteriorSamples[ii, 'FOI'], n)
    ldata$lifespanR <- rep(posteriorSamples[ii,'lifespanR'], n)
    ldata$pref <- rep(posteriorSamples[ii,'bitingPref'], n)
    ldata$disIndMortalityR <- 0
    
  } else if(scenario.parameters == "S1") {
    
    abundanceReservoir <- "blackbird"
    introduction <- "national"
    
    ldata$vertical_transm <- rep(posteriorSamples[ii, 'vertical_transm'], n)
    ldata$transmissionProbHM <- rep(posteriorSamples[ii, 'transmissionProbHM'], n)
    scalingParameter <- posteriorSamples[ii,'scalingParameter']
    ldata$disIndMortalityH <- rep(posteriorSamples[ii,'disIndMortality'], n)
    FOI <- rep(posteriorSamples[ii, 'FOI'], n)
    ldata$lifespanR <- rep(posteriorSamples[ii,'lifespanR'], n)
    ldata$pref <- rep(posteriorSamples[ii,'bitingPref'], n)
    ldata$disIndMortalityR <- 0
    
  } else if(scenario.parameters == "S2") {
    abundanceReservoir <- "uniform"
    introduction <- "south"
    
    ldata$vertical_transm <- rep(posteriorSamples[ii, 'vertical_transm'], n)
    ldata$transmissionProbHM <- rep(posteriorSamples[ii, 'transmissionProbHM'], n)
    scalingParameter <- posteriorSamples[ii,'scalingParameter']
    ldata$disIndMortalityH <- rep(posteriorSamples[ii,'disIndMortality'], n)
    FOI <- rep(posteriorSamples[ii, 'FOI'], n)
    ldata$lifespanR <- rep(posteriorSamples[ii,'lifespanR'], n)
    ldata$pref <- rep(posteriorSamples[ii,'bitingPref'], n)
    ldata$disIndMortalityR <- 0
    
  } else if (scenario.parameters == "S3") {
    abundanceReservoir <- "blackbird"
    introduction <- "south"
    
    ldata$vertical_transm <- rep(posteriorSamples[ii, 'vertical_transm'], n)
    ldata$transmissionProbHM <- rep(posteriorSamples[ii, 'transmissionProbHM'], n)
    scalingParameter <- posteriorSamples[ii,'scalingParameter']
    ldata$disIndMortalityH <- rep(posteriorSamples[ii,'disIndMortality'], n)
    FOI <- rep(posteriorSamples[ii, 'FOI'], n)
    ldata$lifespanR <- rep(posteriorSamples[ii,'lifespanR'], n)
    ldata$pref <- rep(posteriorSamples[ii,'bitingPref'], n)
    ldata$disIndMortalityR <- 0
    
  } else {
    
    print("invalid scenario parameter")
    stop()
    
  }
  

# * run model -------------------------------------------------------------
# different per scenario, depending on one or two-host model and depending on number of parameters in posterior
    
  if (scenario.parameters == "A") {
    # * one-host model
    ## remove the events related to reservoir ageing (using select value 10) to be able to use the events in the one-host system
    events[[1]] <- events[[1]][events[[1]]$select < 8, ]
    
    model <- buildModel(transitions = transitions, compartments = compartments, ldata = ldata, E = E, N = N, v0 = v0, tspan = 92:2550, 
                        startingabundance = startingabundance, introduction = introduction, scalingParameter = scalingParameter, events = events, 
                        pts_fun = pts_fun, age_structure="yes", runtype ="mean", FOI=FOI, useFOI = "yes", verbose=0) 
    
    results <- run(model)
    listofresults[[ii]] <- list(modeloutput = trajectory(results),
                                parameters = posteriorSamples[ii,3:7]) 
    
  } else if (scenario.parameters == "B") {
    # two-host model
    model <- buildModel(transitions = transitions, compartments = compartments, ldata = ldata, E = E, N = N, v0 = v0, tspan = 92:2500, 
                        startingabundance = startingabundance, introduction = introduction, scalingParameter = scalingParameter, events = events, 
                        pts_fun = pts_fun, age_structure="yes", runtype ="mean", FOI=FOI, useFOI = "yes", 
                        abundanceReservoir = abundanceReservoir, verbose=0)    
    
    results <- run(model)
    listofresults[[ii]] <- list(modeloutput = trajectory(results),
                                parameters = posteriorSamples[ii,3:8]) 
    
  } else {
    print(abundanceReservoir)
    
    # two-host model
    model <- buildModel(transitions = transitions, compartments = compartments, ldata = ldata, E = E, N = N, v0 = v0, tspan = 92:2500, 
                        startingabundance = startingabundance, introduction = introduction, scalingParameter = scalingParameter, events = events, 
                        pts_fun = pts_fun, age_structure="yes", runtype ="mean", FOI=FOI, useFOI = "yes", 
                        abundanceReservoir = abundanceReservoir, verbose=0)    
    
    results <- run(model)
    listofresults[[ii]] <- list(modeloutput = trajectory(results),
                                parameters = posteriorSamples[ii,3:9]) 
  }
  
}

save(listofresults, file=filename)

