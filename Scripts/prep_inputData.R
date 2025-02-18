#________________________________________________________________________________
## Title:
# Data Preparation 

## Description:
# This files processes data so it can be used for the analysis/paper called: 

# 1. Create spatial dataframe
# 2. Add columns for mosquito abundance as placeholder
# 3. Estimate bird abundance
# 4. Extract temperature
# 5. Estimate mosquito abundance & emergence per time step

#_________________________________________________________________________________

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
{
  # for bird abundance
  library(raster) 
  library(tidyverse)
  library(tidyr)
  library(terra)
  
  # for temperature
  library(tidync) # package to deal with netcdf files
  library(sf)     # spatial functions
  library(tidyverse)
  library(raster)
  
  # for mosquito
  library(zoo) # to take weekly average
  
  # find day number
  library(lubridate)
  
  # for scales in maps
  library(scales)
}

source("functions_inputData.R")

modelraster <- raster("../Data/model_input/spatialDF5kmNoWater.tif")
modelrasterDF <- read.csv("../Data/model_input/spatialDF5kmNoWater.csv")

# Temperature -------------------------------------------------------------
# Extract temperature from KNMI netcdf files & add to inputData 

# Loop through all files, return dataframe with temperature values (row=raster location, column=date)
# takes a couple of minutes to run

# adds years one by one, need to manually specify the year
filename.temp <- "INTER_OPER_R___TG1_____L3__"
files.temp <- dir("../Data/TemperatureKNMI/2022/")
files.temp <- files.temp[grep(filename.temp, files.temp, TRUE)]

outputfile <- modelrasterDF
setwd("../Data/TemperatureKNMI/2022")

temperature_2022 <- create.temperature(temperaturemaps = files.temp, modelraster = modelraster, outputfile = outputfile)


# add new year to previously created temperature DF
load("C:/Repos/vbd-siminf/Data/model_input/temperature5kmNoWater.RData")
temperature_2022_select <- temperature_2022[c(4:368)] # 369 for leap year, 368 for normal year
temperature <- cbind(temperature, temperature_2022_select)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
save(temperature,file="../Data/model_input/temperature5kmNoWater.RData")


# * plot ------------------------------------------------------------------

ggplot(temperature,aes(x,y, fill=`2019-04-01`)) +
  geom_tile()

# * create plot with mean monthly temperature per year
temperature_forplot <- temperature %>% dplyr::select(-c(x, y, layer))
temperature_forplot <- data.frame(t(temperature_forplot))

temperature_forplot <- temperature_forplot %>%
  rowwise(.) %>%
  mutate(daily.mean.temp = mean(c_across(1:1398))) 

temperature_forplot.date <- temperature_forplot 
temperature_forplot.date$date <- seq(from = as.Date("2016-01-01"), to = as.Date("2022-12-31"), by = 'day')
temperature_forplot.date <- temperature_forplot.date %>%
  mutate(month = format(as.Date(date, format="%Y-%m-%d"),"%m")) %>%
  mutate(year = format(as.Date(date, format="%Y-%m-%d"),"%Y")) %>%
  group_by(month, year) %>%
  summarise(monthly.mean = mean(daily.mean.temp))

monthly.temp.plot <- ggplot(data=temperature_forplot.date) +
  geom_line(aes(x=month, y=monthly.mean, group=year, colour=year), linewidth=1) +
  ylab("Temperature") + xlab("Month") +
  theme_bw() +
  scale_colour_manual(values = brewer.pal(9, "YlOrRd")[3:9]) +
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16), legend.text = element_text(size=16),
        legend.title = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())


monthly.temp.plot
ggsave(monthly.temp.plot, file="../Output/Plots/InputData/meanmonthlytemp.png", width=7, height=5)

# Bird abundance ----------------------------------------------------------
# Create DF with all bird abundance estimates

filename.bird <- "ref_Merel"
files.bird <- dir("../Data/bird_abundance/")
files.bird <- files.bird[grep(filename.bird, files.bird, TRUE)]

outputfile <- modelrasterDF

setwd("../Data/bird_abundance")

birdabundance <- create.birdabundance(abundancemaps = files.bird, modelraster = modelraster, 
                                      outputfile = outputfile, verbose=0)

colnames(birdabundance) <- c("x", "y", "layer", "mean")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
save(birdabundance, file="../Data/model_input/birdabundance5kmNoWater.RData")


# * plot ------------------------------------------------------------------

birdabundance <- rast("../Data/bird_abundance/ref_Merel.asc")
birdabundance <- as.data.frame(birdabundance, xy=TRUE)

ggplot(birdabundance,aes(x,y, fill=ref_Merel)) +
  geom_tile() +
  scale_fill_gradientn(
    "abundance",
    colours = c("white", "lightblue", "darkolivegreen3", "red", "darkred"),
    values = scales::rescale(c(0, 1, 2, 5, 11)),
    limits = c(0, 11)
  ) +
  ylab("") + xlab("") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        strip.text.x=element_text(size=16), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank())

ggsave(filename=paste0("../Output/Plots/inputData/birdabundance_raw.png", sep=""),
       width=5, height=5, bg="white")

# Mosquito abundance ----------------------------------------------------------
# Create DF with all mosquito abundance estimates
# 1 file for each bootstrap, this file contains all timesteps
# loop through all bootstrap folders
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

year <- 2021

for (year in year) {
  print(year)
  
  folders.mosquito <- dir(paste0("../Data/mosquito_abundance_", year))
  filename.mosquito <- "predictedr_count"
  
  outputfile <- modelrasterDF
  mosquitoabundance <- list()
  
  setwd(paste0("../Data/mosquito_abundance_", year))
  
  for (i in 1:1) { # extract only mean estimates
      
    files.mosquito <- dir(paste0(folders.mosquito[i], sep=""))
    
    setwd(paste0(folders.mosquito[i]))
    
    mosquitoabundance.scenario <- create.mosquitoabundance(abundancemaps = files.mosquito, modelraster = modelraster, 
                                                           outputfile = outputfile, verbose=0)
    
    mosquitoabundance[[i]] <- mosquitoabundance.scenario 
    names(mosquitoabundance)[1] <- "mean" 
    
    setwd(paste0("../../mosquito_abundance_", year))
    
  }
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
save(mosquitoabundance, file = paste0("../Data/model_input/mosquitoabundance5kmNoWater_", year, ".RData"))

# Mosquito death and emergence --------------------------------------------

timespan=92:259 # 1 april to 15 sept (until 30 sept (day 274) is the duration of mosquito abundance data, but diapause starts earlier)

if (year == 2017) {
  timespan <- timespan + 365 
} else if (year == 2018) {
  timespan <- timespan + 365*2
} else if (year == 2019) {
  timespan <- timespan + 365*3
} else if (year == 2020) {
  timespan <- timespan + 365*4+1 # +1 bc 2020 is leap year
} else if (year == 2021) {
  timespan <- timespan + 365*5+1
} else if (year == 2022) {
  timespan <- timespan + 365*6+1
}

load("../Data/model_input/temperature5kmNoWater.RData")
temperature <- temperature[,-c(1:3)]    # take out x, y, layer columns
temperature <- temperature[, timespan]  

DeathRate <- sapply(temperature, DeathRateCalc)


load(paste0("../Data/model_input/mosquitoabundance5kmNoWater_", year, ".RData"))

mosquitoabundance_14d <- list()

for (i in 1:length(mosquitoabundance)) {
  
  mosquitoabundance <- mosquitoabundance[[i]] # select file

  # * Derive absolute emergence ---------------------------------------------
  
  mosquitoabundance_14d.scenario <- create.mosquito_14d(mosquitoabundance=mosquitoabundance, timespan=timespan, 
                                               verbose=0)
  
  mosquitoabundance_14d[[i]] <- mosquitoabundance_14d.scenario
}

names(mosquitoabundance_14d)[1] <- "mean" # first list is mean, all others are bootstraps
save(mosquitoabundance_14d, file=paste0("../Data/model_input/mosquitoabundance_14d_5kmNoWater_", year, ".RData"))

mosquitoesEmerging <- list()
for (i in 1:length(mosquitoabundance_14d)) {

  mosquitoabundance_14d <- mosquitoabundance_14d[[i]] # select file
  
  mosquitoesEmerging.scenario <- create.mosquitoEmergence(mosquitoabundance_14d = mosquitoabundance_14d, 
                                                 timespan=timespan, verbose=0)

  mosquitoesEmerging[[i]] <- mosquitoesEmerging.scenario
  
}
names(mosquitoesEmerging)[1] <- "mean" # first list is mean, all others are bootstraps
save(mosquitoesEmerging, file=paste0("../Data/model_input/mosquitoesEmerging_5kmNoWater_", year, ".RData"))


#* Derive absolute death --------------------------------------------------

load(paste0("../Data/model_input/mosquitoabundance_14d_5kmNoWater_", year, ".RData"))

mosquitoesDying <- list()

for (i in 1:length(mosquitoabundance_14d)) {
  
  mosquitoabundance_14d <- mosquitoabundance_14d[[i]] # select file
  
  mosquitoesDying.scenario <- create.mosquitoDeath(mosquitoabundance_14d = mosquitoabundance_14d, DeathRate=DeathRate)
  
  mosquitoesDying[[i]] <- mosquitoesDying.scenario
}

names(mosquitoesDying)[1] <- "mean" # first list is mean, all others are bootstraps
save(mosquitoesDying, file=paste0("../Data/model_input/mosquitoesDying_5kmNoWater_", year, ".RData"))

matplot(t(mosquitoesEmerging[[1]]),type="l")
matplot(t(mosquitoesDying[[1]]),type="l")


# * plot monthly abundance ------------------------------------------------------------------

# * create plot with mean abundance per year
year = 2016
load(paste0("../Data/model_input/mosquitoabundance5kmNoWater_", year, ".RData"))
ab2016 <- mosquitoabundance[["mean"]]

ab2016 <- ab2016 %>% dplyr::select(-c(x, y, layer))
ab2016 <- data.frame(t(ab2016))
ab2016 <- ab2016 %>%
  rowwise(.) %>%
  mutate(daily.mean.ab = mean(c_across(1:1398))) 

plot(ab2016$daily.mean.ab)

ab2016$date <- seq(from = as.Date("2016-04-01"), to = as.Date("2016-09-15"), by = 'day')
ab2016 <- ab2016 %>%
  mutate(month = format(as.Date(date, format="%Y-%m-%d"),"%m")) %>%
  mutate(year = format(as.Date(date, format="%Y-%m-%d"),"%Y")) %>%
  group_by(month, year) %>%
  summarise(monthly.mean = mean(daily.mean.ab)) %>%
  ungroup()


year = 2017
load(paste0("../Data/model_input/mosquitoabundance5kmNoWater_", year, ".RData"))
ab2017 <- mosquitoabundance[["mean"]]

ab2017 <- ab2017 %>% dplyr::select(-c(x, y, layer))
ab2017 <- data.frame(t(ab2017))
ab2017 <- ab2017 %>%
  rowwise(.) %>%
  mutate(daily.mean.ab = mean(c_across(1:1398))) 

plot(ab2017$daily.mean.ab)

ab2017$date <- seq(from = as.Date("2017-04-01"), to = as.Date("2017-09-15"), by = 'day')
ab2017 <- ab2017 %>%
  mutate(month = format(as.Date(date, format="%Y-%m-%d"),"%m")) %>%
  mutate(year = format(as.Date(date, format="%Y-%m-%d"),"%Y")) %>%
  group_by(month, year) %>%
  summarise(monthly.mean = mean(daily.mean.ab)) %>%
  ungroup()


year = 2018
load(paste0("../Data/model_input/mosquitoabundance5kmNoWater_", year, ".RData"))
ab2018 <- mosquitoabundance[["mean"]]

ab2018 <- ab2018 %>% dplyr::select(-c(x, y, layer))
ab2018 <- data.frame(t(ab2018))
ab2018 <- ab2018 %>%
  rowwise(.) %>%
  mutate(daily.mean.ab = mean(c_across(1:1398))) 

ab2018$date <- seq(from = as.Date("2018-04-01"), to = as.Date("2018-09-15"), by = 'day')
ab2018 <- ab2018 %>%
  mutate(month = format(as.Date(date, format="%Y-%m-%d"),"%m")) %>%
  mutate(year = format(as.Date(date, format="%Y-%m-%d"),"%Y")) %>%
  group_by(month, year) %>%
  summarise(monthly.mean = mean(daily.mean.ab)) %>%
  ungroup()


year = 2019
load(paste0("../Data/model_input/mosquitoabundance5kmNoWater_", year, ".RData"))
ab2019 <- mosquitoabundance[["mean"]]

ab2019 <- ab2019 %>% dplyr::select(-c(x, y, layer))
ab2019 <- data.frame(t(ab2019))
ab2019 <- ab2019 %>%
  rowwise(.) %>%
  mutate(daily.mean.ab = mean(c_across(1:1398))) 

ab2019$date <- seq(from = as.Date("2019-04-01"), to = as.Date("2019-09-15"), by = 'day')
ab2019 <- ab2019 %>%
  mutate(month = format(as.Date(date, format="%Y-%m-%d"),"%m")) %>%
  mutate(year = format(as.Date(date, format="%Y-%m-%d"),"%Y")) %>%
  group_by(month, year) %>%
  summarise(monthly.mean = mean(daily.mean.ab)) %>%
  ungroup()


year = 2020
load(paste0("../Data/model_input/mosquitoabundance5kmNoWater_", year, ".RData"))
ab2020 <- mosquitoabundance[["mean"]]

ab2020 <- ab2020 %>% dplyr::select(-c(x, y, layer))
ab2020 <- data.frame(t(ab2020))
ab2020 <- ab2020 %>%
  rowwise(.) %>%
  mutate(daily.mean.ab = mean(c_across(1:1398))) 

ab2020$date <- seq(from = as.Date("2020-04-01"), to = as.Date("2020-09-15"), by = 'day')
ab2020 <- ab2020 %>%
  mutate(month = format(as.Date(date, format="%Y-%m-%d"),"%m")) %>%
  mutate(year = format(as.Date(date, format="%Y-%m-%d"),"%Y")) %>%
  group_by(month, year) %>%
  summarise(monthly.mean = mean(daily.mean.ab)) %>%
  ungroup()


year = 2021
load(paste0("../Data/model_input/mosquitoabundance5kmNoWater_", year, ".RData"))
ab2021 <- mosquitoabundance[["mean"]]

ab2021 <- ab2021 %>% dplyr::select(-c(x, y, layer))
ab2021 <- data.frame(t(ab2021))
ab2021 <- ab2021 %>%
  rowwise(.) %>%
  mutate(daily.mean.ab = mean(c_across(1:1398))) 

ab2021$date <- seq(from = as.Date("2021-04-01"), to = as.Date("2021-09-15"), by = 'day')
ab2021 <- ab2021 %>%
  mutate(month = format(as.Date(date, format="%Y-%m-%d"),"%m")) %>%
  mutate(year = format(as.Date(date, format="%Y-%m-%d"),"%Y")) %>%
  group_by(month, year) %>%
  summarise(monthly.mean = mean(daily.mean.ab)) %>%
  ungroup()


year = 2022
load(paste0("../Data/model_input/mosquitoabundance5kmNoWater_", year, ".RData"))
ab2022 <- mosquitoabundance[["mean"]]

ab2022 <- ab2022 %>% dplyr::select(-c(x, y, layer))
ab2022 <- data.frame(t(ab2022))
ab2022 <- ab2022 %>%
  rowwise(.) %>%
  mutate(daily.mean.ab = mean(c_across(1:1398))) 

ab2022$date <- seq(from = as.Date("2022-04-01"), to = as.Date("2022-09-15"), by = 'day')
ab2022 <- ab2022 %>%
  mutate(month = format(as.Date(date, format="%Y-%m-%d"),"%m")) %>%
  mutate(year = format(as.Date(date, format="%Y-%m-%d"),"%Y")) %>%
  group_by(month, year) %>%
  summarise(monthly.mean = mean(daily.mean.ab)) %>%
  ungroup()


ab2016_22 <- rbind(ab2016, ab2017, ab2018, ab2019, ab2020, ab2021, ab2022)

monthly.abundance.plot <- 
  ggplot(data=ab2016_22) +
  geom_line(aes(x=month, y=monthly.mean, group=year, colour=year), linewidth=1) +
  ggtitle("Mean monthly mosquito abundance") +
  ylab("Abundance (arbitrary unit)") + xlab("Month") +
  theme_bw() +
  scale_colour_manual(values = brewer.pal(9, "YlOrRd")[3:9]) +
  theme(axis.text=element_text(size=16), axis.title = element_text(size=16), legend.text = element_text(size=16),
        legend.title = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  

ggsave(monthly.abundance.plot, file="../Output/Plots/InputData/mosquito/meanmonthlyabundance.png", width=7, height=5)


                                             