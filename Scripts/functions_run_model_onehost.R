#________________________________________________________________________________

## Title:
# functions_main_simulation

## Description:
# This files runs the simulation model for the analysis/paper called: 

#________________________________________________________________________________


buildModel <- function(transitions, compartments, ldata, E, N, v0, tspan, startingabundance, introduction,
                       scalingParameter, events, pts_fun, age_structure = c("yes", "no"), runtype = c("mean", "uncertainty"), 
                       FOI, useFOI = c("yes", "no"), verbose=0) {
  
  if (verbose==1) browser()
  # FOI ---------------------------------------------------------------------

  if (useFOI == 'yes') {
    # use estimated FOI to set prop Rh
    propRh <- 1-exp(-FOI) 
  } else { 
    propRh <- 0
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


    # * initial infections ----------------------------------------------------
    
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
    

    # * blackbird birth ---------------------------------------------------------
    
    # Create bird birth, because dependent on scaling parameter
    birthRate <- 2.5 # 2.5 juveniles per year per adult
    nr.births <- startingabundance[['mean']]$blackbird_abundance * scalingParameter * birthRate
    
    temp <- data.frame(matrix(nrow = 0, ncol = 8)) 
    years <- c(2016, 2017, 2018, 2019, 2020, 2021, 2022)
    
    for (year in years) {

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
    

    # * combine events --------------------------------------------------------

    events <- rbind(events, eventsBirdBirth)
    rm(temp); rm(seeding); rm(eventsBirdBirth)    

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
                       DIha = round(1*startingprevalence$proportionDh * startingabundance[['mean']]['blackbird_abundance'] * scalingParameter))
      
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

