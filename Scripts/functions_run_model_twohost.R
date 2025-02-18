#________________________________________________________________________________

## Title:
# functions_main_simulation

## Description:
# This files runs the simulation model for the analysis/paper called: 

#________________________________________________________________________________


buildModel <- function(transitions, compartments, ldata, E, N, v0, tspan, startingabundance, introduction,
                       scalingParameter, events, pts_fun, age_structure = c("yes", "no"), 
                       runtype = c("mean", "uncertainty"), 
                       FOI, useFOI = c("yes", "no"), 
                       abundanceReservoir = c("uniform", "blackbird"), verbose=0) {
  
  if (verbose==1) browser()
  # FOI ---------------------------------------------------------------------

  if (useFOI == 'yes') {
    # use estimated FOI to set prop Rh
    propRh <- 1-exp(-FOI) 
  } else { 
    propRh <- 0
  }
  

  # Starting abundance reservoir --------------------------------------------
  if (abundanceReservoir == 'blackbird'){
    startingabundance[['mean']]['reservoir_abundance'] <- startingabundance[['mean']]['blackbird_abundance']
  } else if (abundanceReservoir == 'uniform') {
    startingabundance[['mean']]['reservoir_abundance'] <- mean(startingabundance[['mean']]['blackbird_abundance']$blackbird_abundance)
  }
  
  # Starting prevalence -----------------------------------------------------

  startingprevalence <- data.frame(
    proportionSm = 1, 
    proportionEm = 0,
    proportionIm = 0, 
    proportionSh = 1-propRh,
    proportionEh = 0,
    proportionIh = 0,
    proportionRh = propRh,
    proportionDh = 0)
  

  if (runtype == 'mean') {
    
    # Events ------------------------------------------------------------------

    # Select events data frame
    events = events[[1]] # call 'mean'


    # * Initial infections ---------------------------------------------------------------

    if (introduction == "south") {
      
      # create seeding events (different each run)
      seeding <- data.frame(event      = "intTrans",
                            time       = sample(x=seq(from=92, to=290, by=1), size=66), # once every 3 days (when only in South)
                            node       = round(runif(n=66, min=nrow(ldata)*(2/3), max=nrow(ldata)), digits=0), # only in South
                            dest       = 0,
                            n          = 1,
                            proportion = 0,
                            select     = 1,
                            shift      = 3)
      
    } else if (introduction == "national") {
      
      # create seeding events (different each run)
      seeding <- data.frame(event      = "intTrans",
                            time       = seq(from=92, to=290, by=1), # each day
                            node       = round(runif(n=199, min=1, max=nrow(ldata)), digits=0), # random locations
                            dest       = 0,
                            n          = 1,
                            proportion = 0,
                            select     = 1,
                            shift      = 3)
    } else {
      print("no valid introduction scenario")
      stop()
    }
    
    events <- rbind(events, seeding)

    # * Blackbird birth ---------------------------------------------------------

    # Create blackbird birth here, because dependent on scaling parameter & should be fixed number of births
    birthRate <- 2.5 # 2.5 juveniles per year per adult (based in Leslie matrix calculation using juvenile and adult survival)
    nr.births <- startingabundance[['mean']]$blackbird_abundance * scalingParameter * birthRate
    
    temp <- data.frame(matrix(nrow = 0, ncol = 8)) 
    years <- c(2016, 2017, 2018, 2019, 2020, 2021, 2022)
    
    for (year in years) {
      # print(year)
      
      nodes <- length(nr.births)
      timespan <- 122:167 # 1 May to 15 June
      
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
      
      eventsBirdBirth <- data.frame(event      = rep("enter", nodes*length(timespan)),
                                    time       = rep(timespan, each = nodes),
                                    node       = rep(1:nodes, length(timespan)),
                                    dest       = 0,
                                    n          = round(nr.births/length(timespan), digits=0), 
                                    proportion = 0,
                                    select     = 8,
                                    shift      = 0)
      
      temp <- rbind(temp, eventsBirdBirth)
    }
    
    eventsBirdBirth <- temp
    

    # * Reservoir births ------------------------------------------------------

    # Create reservoir birth here, because dependent on scaling parameter & should be fixed number of births.
    # First need to calculate birth rate from Leslie matrix:
      # juvenile mortality rate same as blackbirds
    juvenileSurvival.yearly <- 0.13
      # adult mortality rate based on lifespanR
    adultMortality.daily <- (1/ldata$lifespanR)[1]
    adultSurvival.yearly <- (1-adultMortality.daily)^365
      # use optimise to find the birth rate that leads to an eigenvalue of 1
    calculate_birthRes <- function(x) {
      matrix <- as.matrix(data.frame(c(0, juvenileSurvival.yearly), c(x, adultSurvival.yearly)))
      
      # Calculate the eigenvalues
      eigenvalues <- eigen(matrix)$values
      
      # Calculate the difference from 1
      difference <- eigenvalues - 1
      
      # Return the minimum absolute difference from 1
      return(min(abs(difference)))
    }
    
    # Use optimize to find the value of x that makes an eigenvalue 1
    solution <- optimize(calculate_birthRes, c(0, 10))  
    birthRate.res <- solution$minimum
    print(paste0("birth rate res", birthRate.res, sep=":"))
    

    # add to events, similar to blackbird
    nr.births.res <- startingabundance[['mean']]$reservoir_abundance * scalingParameter * birthRate.res
    
    temp <- data.frame(matrix(nrow = 0, ncol = 8)) 
    years <- c(2016, 2017, 2018, 2019, 2020, 2021, 2022)
    
    for (year in years) {
      # print(year)
      
      nodes <- length(nr.births.res)
      timespan <- 122:167 # 1 May to 15 June
      
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
      
      eventsReservoirBirth <- data.frame(event      = rep("enter", nodes*length(timespan)),
                                         time       = rep(timespan, each = nodes),
                                         node       = rep(1:nodes, length(timespan)),
                                         dest       = 0,
                                         n          = round(nr.births.res/length(timespan), digits=0), 
                                         proportion = 0,
                                         select     = 9, 
                                         shift      = 0)
      
      temp <- rbind(temp, eventsReservoirBirth)
    }
    
    eventsReservoirBirth <- temp
    
    # * combine events --------------------------------------------------------
    
    events <- rbind(events, eventsBirdBirth, eventsReservoirBirth)
    rm(temp); rm(seeding); rm(eventsBirdBirth); rm(eventsReservoirBirth)    

    # Starting conditions -----------------------------------------------------

    # Create starting conditions (u0)
    if (age_structure == 'no') {
    u0 <- data.frame(Sm = round(startingprevalence$proportionSm * startingabundance[['mean']]['mosquito_abundance']),
                     Em = round(startingprevalence$proportionEm * startingabundance[['mean']]['mosquito_abundance']),
                     Im = round(startingprevalence$proportionIm * startingabundance[['mean']]['mosquito_abundance']),
                     Sh = round(startingprevalence$proportionSh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                     Eh = round(startingprevalence$proportionEh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                     Ih = round(startingprevalence$proportionIh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                     Rh = round(startingprevalence$proportionRh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                     Dh = round(startingprevalence$proportionDh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter))
    
    } else if (age_structure == 'yes') {
      
      u0 <- data.frame(Sm = round(startingprevalence$proportionSm * startingabundance[['mean']]['mosquito_abundance']),
                       Em = round(startingprevalence$proportionEm * startingabundance[['mean']]['mosquito_abundance']),
                       Im = round(startingprevalence$proportionIm * startingabundance[['mean']]['mosquito_abundance']),
                       Shj = round(0*startingprevalence$proportionSh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       Ehj = round(0*startingprevalence$proportionEh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       Ihj = round(0*startingprevalence$proportionIh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       Rhj = round(0*startingprevalence$proportionRh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       Dhj = round(0*startingprevalence$proportionDh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       Ddhj = round(0*startingprevalence$proportionDh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       DIhj = round(0*startingprevalence$proportionDh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       Sha = round(1*startingprevalence$proportionSh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       Eha = round(1*startingprevalence$proportionEh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       Iha = round(1*startingprevalence$proportionIh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       Rha = round(1*startingprevalence$proportionRh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       Dha = round(1*startingprevalence$proportionDh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       Ddha = round(1*startingprevalence$proportionDh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       DIha = round(1*startingprevalence$proportionDh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter),
                       Srj = round(0*startingprevalence$proportionSh * startingabundance[['mean']]['reservoir_abundance'] * scalingParameter),
                       Erj = round(0*startingprevalence$proportionEh * startingabundance[['mean']]['reservoir_abundance'] * scalingParameter),
                       Irj = round(0*startingprevalence$proportionIh * startingabundance[['mean']]['reservoir_abundance'] * scalingParameter),
                       Rrj = round(0*startingprevalence$proportionRh * startingabundance[['mean']]['reservoir_abundance'] * scalingParameter),
                       Drj = round(0*startingprevalence$proportionDh * startingabundance[['mean']]['reservoir_abundance'] * scalingParameter),
                       Sra = round(1*startingprevalence$proportionSh * startingabundance[['mean']]['reservoir_abundance'] * scalingParameter),
                       Era = round(1*startingprevalence$proportionEh * startingabundance[['mean']]['reservoir_abundance'] * scalingParameter),
                       Ira = round(1*startingprevalence$proportionIh * startingabundance[['mean']]['reservoir_abundance'] * scalingParameter),
                       Rra = round(1*startingprevalence$proportionRh * startingabundance[['mean']]['reservoir_abundance'] * scalingParameter),
                       Dra = round(1*startingprevalence$proportionDh * startingabundance[['mean']]['reservoir_abundance'] * scalingParameter))
      
    }
  } else if (runtype == 'uncertainty') {
    # draw random dataset. start after 'mean'.
    randomNumber <- sample(2:101, size=1, replace=T) 
    u0 = u0[randomNumber] 
    events = events[randomNumber]
  }
  
  names(u0) <- compartments

  # Mparse ------------------------------------------------------------------

  model <- mparse(transitions = transitions,
                  compartments = compartments,
                  ldata = ldata,
                  events = events,
                  E = E,
                  N = N, 
                  u0 = u0,
                  v0 = v0,
                  tspan = tspan,
                  pts_fun = pts_fun) 
  
  return(model)
}
