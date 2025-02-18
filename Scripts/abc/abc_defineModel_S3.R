
# Load files ------------------------------------------------------------------

load("../Data/model_input/startingabundance.RData") # starting abundance

# events dataframes per dispersal scenario
load("../Data/model_input/events_5kmNoWater_2016_22_agestruc_dispersal.bb.RData") # events bb dispersal

# multi-host model
source("prep_parameters_twohost.R")     # loads compartments, transitions, parameters, initial prevalence, ldata and v0
source("functions_run_model_twohost.R")          # loads 'buildModel' function
source("abc_functions_summaryStatistics.R") # loads functions to calculate summary statistics for ABC


# # load summary statistics from data
# load("../Data/abc_summaryStatistics/deadbird_year_multiyear.Rda")   # vector of relative number of dead bird per month & location
# load("../Data/abc_summaryStatistics/PCRprev_likelihood.Rda")        # prevalence in region N/M/S across year
# load("../Data/abc_summaryStatistics/sero_likelihood.Rda")  # observation process sero data
# load("../Data/abc_summaryStatistics/dead.prev_likelihood.Rda")      # prevalence in region N/M/S across year
# load("../Data/abc_summaryStatistics/popsize.Rda")      # prevalence in region N/M/S across year

# load("../Data/abc_summaryStatistics/node_key.Rda")                  # to link index number of 5x5km grid to index number of 25x25km grid


# function that runs model and calculates summary statistics. 
# parameters to estimate needed as argument

calculate.ss_sim <- function(params, verbose=0) {
  
  if (verbose==1) {browser()}
  
  ### 1a. specify which dispersal scenario and therefore which ldata and events version to use
  
  ldata <- ldata_dispersal.sens
  events <- events_dispersal.bb
  rm(ldata_dispersal.bb, ldata_dispersal.sc20, ldata_dispersal.sc40, ldata_dispersal.sc80, ldata_dispersal.sens, events_dispersal.bb)
  
  ### 1b. specify which distribution to assume for the reservoir population: blackbird or uniform
  
  abundanceReservoir <- "blackbird"
  
  ### 1c. specify with introduction scenario to use: national or south
  
  introduction <- "south"
  
  ### 2. specify parameters to estimate
  
  ldata$vertical_transm <- params[['vertical_transm']]
  ldata$transmissionProbHM <- params[['transmissionProbHM']]
  scalingParameter <- params[['scalingParameter']]
  ldata$disIndMortalityH <- params[['disIndMortality']]
  FOI <- params[['FOI']]
  ldata$lifespanR <- params[['lifespanR']]
  ldata$pref <- params[['bitingPref']]
  # ldata$birthRateR <- (1/params[['lifespanR']])/(46/365) 
  
  print(params)
  print(paste0("lifespan res", mean(ldata$lifespanR), sep=":"))
  # print(paste0("birth rate res", mean(ldata$birthRateR), sep=":"))
  
  ### 3. run model
  model <- buildModel(transitions = transitions, compartments = compartments, ldata = ldata, E = E, N = N, v0 = v0, 
                      tspan=92:2557, startingabundance = startingabundance, introduction = introduction, scalingParameter = scalingParameter, 
                      events = events, pts_fun = pts_fun, age_structure = "yes", runtype ="mean", FOI = FOI, useFOI="yes",
                      abundanceReservoir = abundanceReservoir) 
  
  results <- trajectory(run(model))
  
  ### 4 calculate ss_sim from data
  traject.withindex <- prep.modeloutput(modeloutput=results)
  rm(results)
  ss_sim1 <-  calculate.deadbird_year.pattern(traject.withindex=traject.withindex, output.type="statistic")
  ss_sim2 <- calculate.prevalence.likelihood(traject.withindex=traject.withindex, PCRdataset=PCRprev_likelihood, output.type="statistic")
  ss_sim3 <- calculate.seroprev.likelihood(traject.withindex=traject.withindex, serodataset=sero_likelihood, output.type="statistic")
  ss_sim4 <- calculate.deadprevalence.likelihood(traject.withindex=traject.withindex, deaddataset=dead.prev_likelihood, output.type="statistic")
  ss_sim5 <- calculate.birdpop(traject.withindex=traject.withindex, output.type="statistic")
  
  ss_sim <- list(ss_sim1, ss_sim2, ss_sim3, ss_sim4, ss_sim5)
  rm(traject.withindex)
  
  return(ss_sim)
  
}

# function to compute distance between simulated and observed summary statistics
calculate.distance <- function(ss_sim, ss_obs) {
  # # browser()
  dist1 = sum((ss_sim[[1]]-ss_obs[[1]])^2) # dead bird year
  dist2 = -sum(dbinom(x=PCRprev_likelihood$nr.pos, size=PCRprev_likelihood$nr.samples, prob=ss_sim[[2]],log=TRUE)) # likelihood prevalence
  dist3 = -sum(dbinom(x=sero_likelihood$nr.pos, size=sero_likelihood$nr.samples, prob=ss_sim[[3]],log=TRUE)) # lik seroprev
  dist4 = -sum(dbinom(x=dead.prev_likelihood$nr.pos, size=dead.prev_likelihood$nr.samples, prob=ss_sim[[4]],log=TRUE)) # likelihood prevalence dead
  dist5 = sum((ss_sim[[5]]-ss_obs[[5]])^2) # host popsize 
  
  dist = c(dist1, dist2, dist3, dist4, dist5)
  # print(paste0("distance", seq(1:length(dist)), ":", dist))
  
  return(dist)
}

# combine calculation of summary statistics and of distance to get function for abc input

abc.model <- function(params, ss_obs) {
  ss_sim <- calculate.ss_sim(params=params)
  # print(head(ss_sim))
  dist <- calculate.distance(ss_sim=ss_sim, ss_obs=ss_obs)
  print(paste0("distance", seq(1:length(dist)), ":", dist))
  return(dist)
}



# * model comparison ------------------------------------------------------
MODEL_LIST = list("baseline" = abc.model)


