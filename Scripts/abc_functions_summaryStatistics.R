#________________________________________________________________________________

## Title:
# functions_summary_statistics

## Description:
# This files contains functions to calculate summary statistics from model runs 

#________________________________________________________________________________


prep.modeloutput <- function(modeloutput, verbose=0) {
  
  setDT(modeloutput)
  traject.withindex <- modeloutput[, date := as.Date(time, origin = "2015-12-31")
  ][, year := lubridate::year(date)
  ][, month := lubridate::month(date)
  ][, region := case_when(node <= 466  ~ 1,
                          node > 466 & node <= 932 ~ 2,
                          node > 932 ~ 3)]
    
  return(traject.withindex)
    
}
# Dead bird between-year --------------------------------------------------

# Calculate the number of dead birds per location per year (aggregated locations / larger grid sizes)
# Normalise number per year by total across the whole study period for each location
# This way we use the temporal trend without using the absolute numbers and also capture spatial differences between years

calculate.deadbird_year.pattern <- function(traject.withindex, output.type, verbose=0) {
  if (verbose==1) browser()
  
  # select variables, add year variable
  model.results.withindex <- traject.withindex[,.(node, region, date, year, Duhj, Ddhj, DIhj, Duha, Ddha, DIha)
  ][, Dh.tot := Ddhj + DIhj + Ddha + DIha
  ][!is.na(region),]
  
  # count relative number of birds by year. only include April to October
  model.results.withindex <- model.results.withindex[date == "2016-04-01" | date == "2016-10-31" |
                                                       date == "2017-04-01" | date == "2017-10-31" |
                                                       date == "2018-04-01" | date == "2018-10-31" |
                                                       date == "2019-04-01" | date == "2019-10-31" |
                                                       date == "2020-04-01" | date == "2020-10-31" |
                                                       date == "2021-04-01" | date == "2021-10-31" |
                                                       date == "2022-04-01" | date == "2022-10-31", ] 
  
  dataset.aggregated <- rollup(model.results.withindex, j=sum(Dh.tot), by=c("region", "year", "date"))
  dataset.aggregated <- dataset.aggregated[complete.cases(dataset.aggregated)]
  
  # count number of dead birds that died in each year (so difference between Oct and April of each year)
  dataset.aggregated <- dataset.aggregated %>%
    rename(Dh.tot = V1) %>%
    group_by(region) %>% 
    mutate(Dh_year := Dh.tot-lag(Dh.tot)) %>%
    filter(date == "2016-10-31" | date == "2017-10-31" | date == "2018-10-31" | date == "2019-10-31" | 
             date == "2020-10-31" | date == "2021-10-31" | date == "2022-10-31") %>%
    ungroup()
  
  # relative number of birds per cell
  dataset.aggregated <- dataset.aggregated %>%
    group_by(region) %>%
    mutate(total.index = sum(Dh_year)) %>%
    ungroup %>%
    rowwise() %>%
    mutate(relative.deadBirds = Dh_year/total.index) %>%
    arrange(region)
  
  if (output.type == "statistic") {
    setDT(dataset.aggregated)
    statistic <- dataset.aggregated[,.(relative.deadBirds)]
    
    statistic <- as.vector(statistic[[1]])
    
    return(statistic)
  }
  
  if (output.type == "table") {
    dataset.aggregated <- dataset.aggregated %>%
      dplyr::select(region, year, relative.deadBirds)
    
    return(dataset.aggregated)
  }
  
  if (output.type == "plot") {
    plot.deadbird_region <-
      ggplot(data=dataset.aggregated, aes(x=factor(region), y=relative.deadBirds, fill=factor(year))) +
      geom_col(position=position_stack(reverse = TRUE)) +
      scale_fill_viridis_d() +
      ylab("Proportion of reported dead birds") + labs(fill="Year") +
      ggtitle("Dead birds (untested)")   +
      scale_x_discrete(name= "Region", labels=c("North", "Middle", "South"))  +
      theme_bw() + 
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    return(plot.deadbird_region)
  }
}


# Live prevalence ---------------------------------------------

calculate.prevalence.likelihood <- function(traject.withindex, PCRdataset, output.type, verbose=0) {
  if (verbose==1) browser()
  
  # select variables
  traject.withindex <- traject.withindex[!is.na(region),
  ][,.(Shj, Ehj, Ihj, Rhj, Sha, Eha, Iha, Rha, region, date, year, month)]
  
  # aggregate to region 
  traject.largegrid <- rollup(traject.withindex, j=lapply(.SD,sum), by=c("region", "date", "year", "month"))
  traject.largegrid <- traject.largegrid[complete.cases(traject.largegrid)]
  
  # aggregate to month
  traject.largegrid <- rollup(traject.largegrid, j=lapply(.SD,mean), by=c("region", "year", "month"))
  
  # calculate prevalence each month and region
  output.month <- traject.largegrid[, monthmean := (Ehj+Eha+Ihj+Iha)/(Shj+Ehj+Ihj+Rhj+Sha+Eha+Iha+Rha)]
  
  ll.input <- output.month %>%
    arrange(month, year, region) 
  
  # only keep months with observations
  ll.input.reduced <- inner_join(ll.input, PCRdataset, by=c("month", "year", "region"))
  
  if (output.type == "statistic") {
    ll.input.reduced <- ll.input.reduced %>%
      dplyr::select(monthmean)
    
    # replace prevalence = 0 with prevalence = 1/10000000 to avoid loglik of -Inf
    ll.input.reduced <- ll.input.reduced %>%
      mutate(monthmean = replace(monthmean, monthmean == 0, 1/10000000)) %>%
      mutate(monthmean = replace(monthmean, is.nan(monthmean), 1/10000000))
    
    statistic <- as.vector(ll.input.reduced[[1]]) 
    
    return(statistic)
  }
  
  if (output.type == "table") {
    plot.year <- ll.input.reduced %>% group_by(region, year) %>%
      summarise(yearmean = mean(monthmean))
    return(plot.year)
  }
  
  if (output.type == "plot") {
    plot.year <- ll.input.reduced %>% group_by(region, year) %>%
      summarise(yearmean = mean(monthmean))
    
    plot.liveprev.region <- 
      ggplot(data=plot.year, aes(x=as.factor(region), y=yearmean, fill=factor(year))) +
      facet_wrap(vars(year), nrow=3) +
      geom_col() +
      ylim(0, 0.5) + ylab("Prevalence") +
      scale_x_discrete(name= "Region", labels=c("North", "Middle", "South")) +
      ggtitle("Prevalence live birds") +
      scale_fill_viridis_d() +
      theme_bw() + 
      theme(legend.position = "none", panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    return(plot.liveprev.region)
  }
  
}


# Seroprev ----------------------------------------------------------

calculate.seroprev.likelihood <- function(traject.withindex, serodataset, output.type, verbose=0) {
  if (verbose==1) browser()
  
  # select variables
  traject.withindex <- traject.withindex[!is.na(region),
  ][,.(Shj, Ehj, Ihj, Rhj, Sha, Eha, Iha, Rha, region, date, year, month)]
  
  # aggregate to region 
  traject.largegrid <- rollup(traject.withindex, j=lapply(.SD,sum), by=c("region", "date", "year", "month"))
  traject.largegrid <- traject.largegrid[complete.cases(traject.largegrid)]
  
  # aggregate to month
  traject.largegrid <- rollup(traject.largegrid, j=lapply(.SD,mean), by=c("region", "year", "month"))
  
  # calculate seroprevalence each month and region
  output.month <- traject.largegrid[, monthmean := (Rhj+Rha)/(Shj+Ehj+Ihj+Rhj+Sha+Eha+Iha+Rha)]
  
  ll.input <- output.month %>%
    arrange(month, year, region) 
  
  # only keep months with observations
  ll.input.reduced <- inner_join(ll.input, serodataset, by=c("month", "year", "region"))
  
  
  if (output.type == "statistic") {
    ll.input.reduced <- ll.input.reduced %>%
      dplyr::select(monthmean)
    
    # replace seroprevalence = 0 with seroprevalence = 1/10000000 to avoid loglik of -Inf
    ll.input.reduced <- ll.input.reduced %>%
      mutate(monthmean = replace(monthmean, monthmean == 0, 1/10000000)) %>%
      mutate(monthmean = replace(monthmean, is.nan(monthmean), 1/10000000))
    
    statistic <- as.vector(ll.input.reduced[[1]]) 
    
    return(statistic)
  }
  
  if (output.type == "table") {
    plot.year <- ll.input.reduced %>% group_by(region, year) %>%
      summarise(yearmean = mean(monthmean))
    return(plot.year)
  }
  
  if (output.type == "plot") {
    plot.year <- ll.input.reduced %>% group_by(region, year) %>%
      summarise(yearmean = mean(monthmean))
    
    plot.sero.region <- 
      ggplot(data=plot.year, aes(x=as.factor(region), y=yearmean, fill=factor(year))) +
      facet_wrap(vars(year), nrow=3) +
      geom_col() +
      ylim(0,1) + ylab("Seroprevalence") +
      scale_x_discrete(name= "Region", labels=c("North", "Middle", "South")) +
      ggtitle("Seroprevalence live birds") +
      scale_fill_viridis_d() +
      theme_bw() + 
      theme(legend.position = "none", panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    return(plot.sero.region)
  }
  
}

# Dead prevalence ---------------------------------------------

calculate.deadprevalence.likelihood <- function(traject.withindex, deaddataset, output.type, verbose=0) {
  
  # select variables
  traject.withindex <- traject.withindex[!is.na(region),
  ][,.(Ddhj, DIhj, Ddha, DIha, region, date, year, month)]
  
  # aggregate to region 
  traject.largegrid <- rollup(traject.withindex, j=lapply(.SD,sum), by=c("region", "date", "year", "month"))
  traject.largegrid <- traject.largegrid[complete.cases(traject.largegrid)]
  
  if (verbose==1) browser()
  
  # only birds died in the past day can be tested, so only count newly died birds
  traject.largegrid <- traject.largegrid %>%
    group_by(region) %>%
    arrange(date) %>%
    mutate(Ddhj_new = Ddhj - lag(Ddhj, n=1)) %>%
    mutate(DIhj_new = DIhj - lag(DIhj, n=1)) %>%
    mutate(Ddha_new = Ddha - lag(Ddha, n=1)) %>%
    mutate(DIha_new = DIha - lag(DIha, n=1))
  
  # aggregate to month
  setDT(traject.largegrid)
  traject.largegrid <- rollup(traject.largegrid, j=lapply(.SD,mean), by=c("region", "year", "month"))
  
  # calculate prevalence each month and region
  output.month <- traject.largegrid[, monthmean := (DIhj_new+DIha_new)/(Ddhj_new+DIhj_new+Ddha_new+DIha_new)]
  
  ll.input <- output.month %>%
    arrange(month, year, region) 
  
  # only keep weeks/months with observations

  ll.input.reduced <- inner_join(ll.input, deaddataset, by=c("month", "year", "region"))
  
  if (output.type == "statistic") {
    ll.input.reduced <- ll.input.reduced %>%
      dplyr::select(monthmean)
    
    # replace prevalence = 0 with prevalence = 1/10000000 to avoid loglik of -Inf
    ll.input.reduced <- ll.input.reduced %>%
      mutate(monthmean = replace(monthmean, monthmean == 0, 1/10000000)) %>%
      mutate(monthmean = replace(monthmean, is.nan(monthmean), 1/10000000))
    
    statistic <- as.vector(ll.input.reduced[[1]]) 
    
    return(statistic)
  }
  
  if (output.type == "table") {
    plot.year <- ll.input.reduced %>% group_by(region, year) %>%
      summarise(yearmean = mean(monthmean))
    return(plot.year)
  }
  
  if (output.type == "plot") {
    plot.year <- ll.input.reduced %>% group_by(region, year) %>%
      summarise(yearmean = mean(monthmean))
    
    plot.deadprev.region <-
      ggplot(data=plot.year, aes(x=as.factor(region), y=yearmean, fill=factor(year))) +
      facet_wrap(vars(year), nrow=3) +
      geom_col() +
      ylim(0, 1.2) + ylab("Prevalence") +
      scale_x_discrete(name= "Region", labels=c("North", "Middle", "South")) +
      ggtitle("Prevalence dead birds") +
      scale_fill_viridis_d() +
      theme_bw() + 
      theme(legend.position = "none", panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    return(plot.deadprev.region)
  }
  
}

# Bird population size ---------------------------------------------

calculate.birdpop <- function(traject.withindex, output.type, verbose=0) {
  if (verbose==1) browser()
  
  # 1. filter on months in which counts are conducted (March, April, May, June, July)
  traject.withindex <- traject.withindex %>%
    filter(month>3 & month<8)
  
  # 2. select variables
  traject.withindex <- traject.withindex[!is.na(region),
  ][, hostpopsize := Shj+Ehj+Ihj+Rhj+Sha+Eha+Iha+Rha 
  ][,.(hostpopsize, region, year)]
  
  # 3. aggregate to region and year - calculate mean
  traject.largegrid <- rollup(traject.withindex, j=lapply(.SD,mean), by=c("region", "year"))
  traject.largegrid <- traject.largegrid[complete.cases(traject.largegrid)]
  
  r1 <- traject.largegrid %>% filter(region==1) %>% arrange(year) %>%
    mutate(rel.pop = hostpopsize/hostpopsize[1])
  r2 <- traject.largegrid %>% filter(region==2) %>% arrange(year) %>%
    mutate(rel.pop = hostpopsize/hostpopsize[1])
  r3 <- traject.largegrid %>% filter(region==3) %>% arrange(year) %>%
    mutate(rel.pop = hostpopsize/hostpopsize[1])
  
  popsize <- rbind(r1, r2, r3)
  
  if (output.type == "statistic") {
    statistic <- as.vector(popsize$rel.pop) 
    return(statistic)
  }
  
  if (output.type == "table"){
    popsize <- popsize %>% dplyr::select(region, year, rel.pop)
    popsize <- as.data.frame(popsize)
    return(popsize)
  }
  
  if (output.type == "plot") {
    popsize <- popsize %>% dplyr::select(region, year, rel.pop)
    popsize <- as.data.frame(popsize)
    
    plot.birdpop <-
      ggplot(data=popsize, aes(x=year, y=rel.pop, group=as.factor(region), colour=as.factor(region))) +
      geom_line() +
      ylim(0.5, 1.1) + ylab("Relative population size") + xlab("Year") +
      ggtitle("Host population size") +  
      scale_color_brewer(name= "Region", labels=c("North", "Middle", "South"), palette="Dark2") +
      theme_bw() +
      theme(legend.position = c(0.85, 0.8), panel.grid.major = element_blank(),panel.grid.minor = element_blank())
    
    return(plot.birdpop)
  }
}
