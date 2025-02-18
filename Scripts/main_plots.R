# All outputs are created from multiple stochastic runs

# visual fit to data
# R0 over time

rm(list=ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)
library(lubridate)
library(patchwork)
library(scales)


source("functions_processOutput.R") 


# Combine all processed runs -------------------------------------------------------

filename <- "multi_20240805_S1"
runname <- "runs_20240805_S1"

# Number of folders across which processedRuns is spread
num_folders <- 10  
folder_paths <- sprintf(paste0("../Output/Simulations/", filename, "/processedRuns_%d", sep=""), seq(1, num_folders))

# Combine the data from all processedRuns files in the specified folders
ss_sim_posterior <- combine_rds_files(folder_paths = folder_paths, file_name = "summ_stat_posterior.rds")
simOutput <- combine_rds_files(folder_paths = folder_paths, file_name = "simOutput.rds") # select relevant analysis
simOutput <- combine_rds_files(folder_paths = folder_paths, file_name = "simOutput_space.rds")

size <- 13
nr.sim <- 100


# Visual fit patterns ------------------------------------------------------------
dead.untested <- get_ss.sim.output(summ.stat.dimensions = 21, ss_sim_posterior = ss_sim_posterior, patternnr=1)

live.prev <- get_ss.sim.output(summ.stat.dimensions = 21, ss_sim_posterior = ss_sim_posterior, patternnr=2)

sero.prev <- get_ss.sim.output(summ.stat.dimensions = 20, ss_sim_posterior = ss_sim_posterior, patternnr=3)

dead.prev <- get_ss.sim.output(summ.stat.dimensions = 19, ss_sim_posterior = ss_sim_posterior, patternnr=4)

bird.pop <- get_ss.sim.output(summ.stat.dimensions = 21, ss_sim_posterior = ss_sim_posterior, patternnr=5)


# * dead bird untested -------------------------------------------------------
load("../Data/abc_summaryStatistics/forPlotting/data.deadbird.untested.Rda")

dead.untested <- dead.untested %>%   
  pivot_longer(cols=c(1:nr.sim), names_to = "Run", values_to = "Value")

plot.dead.untested <-
  ggplot() +
  geom_violin(data=dead.untested, aes(x=factor(year), y=Value, fill=factor(region)), scale="width") +
  facet_grid(rows=vars(region), labeller = labeller(region = c("1" = "North", "2" = "Middle", "3" = "South")), 
             switch="y") +
  geom_point(data=deadbird_year, aes(x=factor(year), y=relative.reportedBirds, fill=factor(region))) +
  geom_linerange(data=deadbird_year, aes(x=factor(year), ymin=relative.reportedBirds-1.96*se, ymax=relative.reportedBirds+1.96*se)) +
  scale_fill_viridis_d() +
  ylim(0, 0.65) + 
  ylab("Proportion") + 
  xlab("") +
  ggtitle("Relative number of dead blackbirds (untested)") +
  theme_bw() + 
  theme(axis.text=element_text(size=size), 
        axis.title=element_text(size=size+2),
        legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=size), 
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"))

 
  
# * live prevalence  ----------------------------------------------------
load("../Data/abc_summaryStatistics/forPlotting/data.liveprev.Rda")

live.prev <- live.prev %>%   
  pivot_longer(cols=c(1:nr.sim), names_to = "Run", values_to = "Value") 

plot.liveprev.alt <-
  ggplot() +
  geom_violin(data=live.prev, aes(x=as.factor(year), y=Value, fill=factor(region)), scale="width") +
  facet_grid(rows=vars(region), labeller = labeller(region = c("1" = "North", "2" = "Middle", "3" = "South")), 
             switch="y") +
  geom_point(data=prev.region, aes(x=as.factor(year), y=prevalence, fill=factor(region))) +
  geom_linerange(data=prev.region, aes(x=as.factor(year), ymin = ymin, 
                                       ymax = ifelse(prevalence + 1.96*se < 0.3, prevalence + 1.96*se, 0.3))) + 
  ylim(0, 0.3) + 
  ylab("") +
  xlab("") +
    # scale_x_discrete(name= "Year") +
  ggtitle("PCR prevalence live blackbirds") +
  scale_fill_viridis_d() +
  theme_bw() + 
  theme(axis.text=element_text(size=size), 
        axis.title=element_text(size=size+2),
        legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_text(size=size), 
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"))

 
# * sero prev --------------------------------------------------------
load("../Data/abc_summaryStatistics/forPlotting/data.seroprev.Rda")
  
sero.prev <- sero.prev %>%   
  pivot_longer(cols=c(1:nr.sim), names_to = "Run", values_to = "Value")

plot.seroprev.alt <-
ggplot() +
  geom_violin(data=sero.prev, aes(x=as.factor(year), y=Value, fill=factor(region)), scale="width") +
  facet_grid(rows=vars(region), labeller = labeller(region = c("1" = "North", "2" = "Middle", "3" = "South")), 
             switch="y") +
  geom_point(data=sero.region, aes(x=as.factor(year), y=prevalence, fill=factor(region))) +
  geom_linerange(data=sero.region, aes(x=as.factor(year), ymin = ymin, ymax = prevalence + 1.96*se)) +
  ylim(0, 1) + 
  ylab("") +
  xlab("") +
  ggtitle("Seroprevalence live blackbirds") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.text=element_text(size=size), 
        axis.title=element_text(size=size+2),
        legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"))


# * dead prev -----------------------------------------------------
load("../Data/abc_summaryStatistics/forPlotting/data.deadprev.Rda")

dead.prev <- dead.prev %>%   
  pivot_longer(cols=c(1:nr.sim), names_to = "Run", values_to = "Value")

plot.deadprev.alt <-
  ggplot() +
  geom_violin(data=dead.prev, aes(x=as.factor(year), y=Value, fill=factor(region)), scale="width") +
  facet_grid(rows=vars(region), labeller = labeller(region = c("1" = "North", "2" = "Middle", "3" = "South")), 
             switch="y") +
  geom_point(data=dead.prev.region, aes(x=as.factor(year), y=prevalence, fill=factor(region))) +
  geom_linerange(data=dead.prev.region, aes(x=as.factor(year), ymin = ymin, 
                                            ymax = ifelse(prevalence + 1.96*se < 1, prevalence + 1.96*se, 1))) + 
  ylim(0, 1)  + 
  ylab("") +
  xlab("") +
  ggtitle("PCR prevalence dead blackbirds") +
  scale_fill_viridis_d() +
  theme_bw() + 
  theme(axis.text=element_text(size=size), 
        axis.title=element_text(size=size+2),
        legend.position = "none", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_blank(), 
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"))


# * bird pop ----------------------------------------------------------------
load("../Data/abc_summaryStatistics/forPlotting/popsize.Rda")

bird.pop <- bird.pop %>%   
  pivot_longer(cols=c(1:nr.sim), names_to = "Run", values_to = "Value")

plot.pop <- 
ggplot() +
  geom_line(data=popSize, aes(x=as.factor(year), y=pop, group=as.factor(region)), linewidth=0.5) +
  facet_grid(rows=vars(region), labeller = labeller(region = c("1" = "North", "2" = "Middle", "3" = "South")),
             switch="y") +
  geom_violin(data=bird.pop, aes(x=as.factor(year), y=Value, fill=factor(region)), scale="width") +
  geom_point(data=popSize, aes(x=as.factor(year), y=pop, group=as.factor(region)), size=2) +

  ylim(0.5, 1.1) + 
  ylab("Relative population size") + 
  xlab("") +
  ggtitle("Blackbird population size") +  
  scale_color_viridis_d(name= "Region", labels=c("North", "Middle", "South")) +
  scale_fill_viridis_d(name= "Region", labels=c("North", "Middle", "South")) +
  theme_bw() +
  theme(axis.text=element_text(size=size),
        axis.title=element_text(size=size+2),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text=element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"))

# * combine plots ---------------------------------------------------------
fit.plots <- (plot.dead.untested + plot.pop) / (plot.liveprev.alt + plot.seroprev.alt + plot.deadprev.alt) + 
  plot_annotation(tag_levels = 'A') + plot_annotation(tag_levels = 'A')
fit.plots 

ggsave(filename=paste0("../Output/Plots/dataFits/", runname, "/summ_stat_posteriors.png", sep=""), 
       plot=fit.plots, 
       width=14, height=10, bg="white")


# Reproduction number ----------------------------------------------------------------------

## R0
R0.processed <- process.R0(simOutput=simOutput, average="no")
save(R0.processed, file=sprintf(paste("../Output/Simulations/", filename, '/R0.processed.Rdata', sep='')))

## Re
Re.processed <- process.Re(simOutput=simOutput, average="no")
save(Re.processed, file=sprintf(paste("../Output/Simulations/", filename, '/Re.processed.Rdata', sep='')))

# * R0 --------------------------------------------------------------------

R0.processed <- R0.processed %>%
  mutate(year = as.numeric(format(date,'%Y'))) %>%
  mutate(Week = week(date)) %>%
  mutate(doy = as.numeric(strftime(date, format = "%j"))) 

meanR0 <- R0.processed %>% 
  mutate(month = as.numeric(format(date,'%m'))) %>%
  filter(month>=4 & month<=10) %>%
  group_by(year) %>%
  summarise(meanR0 = mean(mean.all))

ggplot() +
  geom_line(data=R0.processed, aes(x=doy, y=mean.all)) + 
  geom_ribbon(data=R0.processed, aes(x=doy , ymin=lb.all, ymax=ub.all), alpha=0.3) +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_text(data=meanR0, aes(x=100, y=13, label=paste("Mean: ", format(round(meanR0, digits=2), nsmall=2), sep="")), size=6) +
  facet_wrap(vars(year), nrow=7) +
  ylab("R0") + xlab("Date") +
  scale_x_continuous(breaks = c(92,153,214,275), 
                     labels = c("April","June","August","October")) +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())

ggsave(filename=paste0("../Output/Plots/dataFits/", runname, "/R0_aprOct.png", sep=""),
width=8, height=11, bg="white")

# ** Maps -----------------------------------------------------------------
R0.processed.space <- process.R0.space(simOutput=simOutput, average="no")
R0.processed.space <- R0.processed.space[1:(7*1398), ] #remove final 1398 rows, because they are the average that should not have been added

save(R0.processed.space, file=sprintf(paste("../Output/Simulations/", filename, '/R0.processed.space.Rdata', sep='')))

# add coordinates
load("../Data/model_input/startingabundance.RData") # starting abundance
spatialInfo <- startingabundance[["mean"]][,1:3]
R0.processed.space.spatial <- R0.processed.space
R0.processed.space.spatial$x <- rep(startingabundance[["mean"]][,1], 7)
R0.processed.space.spatial$y <- rep(startingabundance[["mean"]][,2], 7)

max(R0.processed.space.spatial$mean.all)
min(R0.processed.space.spatial$mean.all)
mean(R0.processed.space.spatial$mean.all)

meanR0 <- R0.processed.space.spatial %>% 
  group_by(year) %>%
  summarise(meanR0 = mean(mean.all))


ggplot(R0.processed.space.spatial,aes(x,y, fill=mean.all))+
  facet_wrap(~year, nrow=2) +
  geom_raster() +
  scale_fill_gradientn(
    "R0",
    colours = c("darkgreen", "lightgreen", "white", "yellow", "red"),
    values = rescale(c(0, 0.99, 1, 3, 8.6)),
    limits = c(0, 8.6)
  ) +
  ylab("") + xlab("") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        strip.text.x=element_text(size=10), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"))

ggsave(paste0(filename="../Output/Plots/dataFits/", runname, "/R0_space.png", sep=""),
       width=6.5, height=4, bg="white")

## aggregate all years
R0.space.allyears <- R0.processed.space.spatial %>% 
  group_by(x,y) %>% summarise(mean.allyears = mean(mean))


# * Partial R0 ------------------------------------------------------------
R0.processed <- R0.processed %>%
  mutate(year = as.numeric(format(date,'%Y'))) %>%
  mutate(Week = week(date)) %>%
  mutate(doy = as.numeric(strftime(date, format = "%j")))

partialR0.year <- R0.processed %>% 
  group_by(year) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols=c(2:1101), names_to = "partialR0", values_to = "Value") %>%
  dplyr::select("year", "partialR0", "Value") %>% 
  mutate(partialR0 = sub("_.*", "", partialR0)) %>%
  filter(partialR0 == "R0.mosq.juv" | partialR0 == "R0.juv.mosq" | partialR0 == "R0.mosq.adu" | 
         partialR0 == "R0.adu.mosq" | partialR0 == "R0.mosq.resjuv" | partialR0 == "R0.resjuv.mosq" | partialR0 == "R0.mosq.resadu" | partialR0 == "R0.resadu.mosq")

partialR0.year <- partialR0.year[partialR0.year$partialR0 != "R0.all", ]
# Convert partialR0 to factor with desired order
partialR0.year$partialR0 <- factor(partialR0.year$partialR0, levels = c("R0.mosq.juv", "R0.juv.mosq", "R0.mosq.adu", "R0.adu.mosq",
                                                                        "R0.mosq.resjuv", "R0.resjuv.mosq", "R0.mosq.resadu", "R0.resadu.mosq"))

partialR0 <-
  ggplot(data=partialR0.year) +
  geom_violin(aes(x=factor(year), y=Value, fill=partialR0), show.legend = FALSE) +
  facet_wrap(~partialR0, scales="free", nrow=4, 
             labeller = labeller(partialR0 = c("R0.adu.mosq" = "Adult blackbird to mosquito", 
                                               "R0.juv.mosq" = "Juvenile blackbird to mosquito", 
                                               "R0.mosq.adu" = "Mosquito to adult blackbird",
                                               "R0.mosq.juv" = "Mosquito to juvenile blackbird",
                                               "R0.mosq.resjuv" = "Mosquito to juvenile reservoir host",
                                               "R0.resjuv.mosq" = "Juvenile reservoir host to mosquito",
                                               "R0.mosq.resadu" = "Mosquito to adult reservoir host",
                                               "R0.resadu.mosq" = "Adult reservoir host to mosquito"))) +
  labs(x = "",
       y = "Partial R0") +
  theme(axis.title.x=element_blank(),
        axis.text.x.bottom = element_blank(),
        legend.position = "none") +
  theme_bw() 


ggsave(filename=paste0("../Output/Plots/dataFits/", runname, "/partialR0_violin.png", sep=""),
       width=6, height=7, bg="white")




# * Re --------------------------------------------------------------------
load(paste0("../Output/Simulations/", filename, "/Re.processed.Rdata", sep=""))

Re.processed <- Re.processed %>%
  mutate(year = as.numeric(format(date,'%Y'))) %>%
  mutate(Week = week(date)) %>%
  mutate(doy = as.numeric(strftime(date, format = "%j")))

meanRe <- Re.processed %>% 
  mutate(month = as.numeric(format(date,'%m'))) %>%
  filter(month>=4 & month<=9) %>%
  group_by(year) %>%
  summarise(meanRe = mean(mean.all))

ggplot() +
  geom_line(data=Re.processed, aes(x=doy, y=mean.all)) + 
  geom_ribbon(data=Re.processed, aes(x=doy , ymin=lb.all, ymax=ub.all), alpha=0.3) +
  geom_hline(yintercept = 1, linetype="dashed") +
  geom_text(data=meanRe, aes(x=100, y=13, label=paste("Mean: ", format(round(meanRe, digits=2), nsmall=2), sep="")), size=6) +
  facet_wrap(vars(year), nrow=7) +
  ylab("Re") + xlab("Date") +
  scale_x_continuous(breaks = c(92,153,214,275), 
                     labels = c("April","June","August","October")) +
  theme_bw() +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), strip.text.x=element_text(size=16),
        legend.title=element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank())


ggsave(filename=paste0("../Output/Plots/dataFits/", runname, "/Re.png", sep=""),
       width=8, height=11, bg="white")

# ** Maps -----------------------------------------------------------------
Re.processed.space <- process.Re.space(simOutput=simOutput, average="no")
Re.processed.space <- Re.processed.space[1:(7*1398), ] #remove final 1398 rows, because they are the average that should not have been added

save(Re.processed.space, file=sprintf(paste("../Output/Simulations/", filename, '/Re.processed.space.Rdata', sep='')))

# add coordinates
load("../Data/model_input/startingabundance.RData") # starting abundance
spatialInfo <- startingabundance[["mean"]][,1:3]
Re.processed.space.spatial <- Re.processed.space
Re.processed.space.spatial$x <- rep(startingabundance[["mean"]][,1], 7)
Re.processed.space.spatial$y <- rep(startingabundance[["mean"]][,2], 7)

max(Re.processed.space.spatial$mean.all)
min(Re.processed.space.spatial$mean.all)
mean(Re.processed.space.spatial$mean.all)

ggplot(Re.processed.space.spatial,aes(x,y, fill=mean.all))+
  facet_wrap(~year, nrow=2) +
  geom_raster() +
  scale_fill_gradientn(
    "Re",
    colours = c("darkgreen", "lightgreen", "white", "yellow", "red"),
    values = rescale(c(0, 0.99, 1, 2, 7.4)),
    limits = c(0, 7.4)
  ) +
  ylab("") + xlab("") +
  theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        strip.text.x=element_text(size=10), 
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"))

ggsave(paste0(filename="../Output/Plots/dataFits/", runname, "/Re_space.png", sep=""),
       width=6.5, height=4, bg="white")

## aggregate all years
Re.space.allyears <- Re.processed.space.spatial %>% 
  group_by(x,y) %>% summarise(mean.allyears = mean(mean))

# * Partial Re ------------------------------------------------------------
Re.processed <- Re.processed %>%
  mutate(year = as.numeric(format(date,'%Y'))) %>%
  mutate(Week = week(date)) %>%
  mutate(doy = as.numeric(strftime(date, format = "%j")))


partialRe.year <- Re.processed %>% 
  group_by(year) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols=c(2:748), names_to = "partialRe", values_to = "Value") %>%
  dplyr::select("year", "partialRe", "Value") %>% 
  mutate(partialRe = sub("_.*", "", partialRe)) %>%
  filter(partialRe == "Re.mosq.juv" | partialRe == "Re.juv.mosq" | partialRe == "Re.mosq.adu" | partialRe == "Re.adu.mosq" | 
           partialRe == "Re.mosq.resjuv" | partialRe == "Re.resjuv.mosq" | partialRe == "Re.mosq.resadu" | partialRe == "Re.resadu.mosq")

partialRe.year <- partialRe.year[partialRe.year$partialRe != "Re.all", ]
# Convert partialRe to factor with desired order
partialRe.year$partialRe <- factor(partialRe.year$partialRe, levels = c("Re.mosq.juv", "Re.juv.mosq", "Re.mosq.adu", "Re.adu.mosq", 
                                                                        "Re.mosq.resjuv", "Re.resjuv.mosq",  "Re.mosq.resadu", "Re.resadu.mosq"))

partialRe <-
  ggplot(data=partialRe.year) +
  geom_violin(aes(x=factor(year), y=Value, fill=partialRe), show.legend = FALSE) +
  facet_wrap(~partialRe, scales="free", nrow=4, 
             labeller = labeller(partialRe = c("Re.adu.mosq" = "Adult blackbird to mosquito", 
                                               "Re.juv.mosq" = "Juvenile blackbird to mosquito", 
                                               "Re.mosq.adu" = "Mosquito to adult blackbird",
                                               "Re.mosq.juv" = "Mosquito to juvenile blackbird",
                                               "Re.mosq.resjuv" = "Mosquito to juvenile reservoir host",
                                               "Re.resjuv.mosq" = "Juvenile reservoir host to mosquito",
                                               "Re.mosq.resadu" = "Mosquito to adult reservoir host",
                                               "Re.resadu.mosq" = "Adult reservoir host to mosquito"))) +
  labs(x = "",
       y = "Partial Re") +
  theme(axis.title.x=element_blank(),
        axis.text.x.bottom = element_blank(),
        legend.position = "none") +
  theme_bw() 

partialRe

ggsave(filename=paste0("../Output/Plots/dataFits/", runname, "/partialRe_violin.png", sep=""),
       width=6, height=7, bg="white")


partialR.plots <- (partialR0 + partialRe) + 
  plot_annotation(tag_levels = 'A') 
partialR.plots 

ggsave(filename=paste0("../Output/Plots/dataFits/", runname, "/partialR_panel.png", sep=""), 
       plot=partialR.plots, 
       width=12, height=8, bg="white")

## violin plots each route Re
partialRe.year <- Re.processed %>% 
  group_by(year) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols=c(2:748), names_to = "partialRe", values_to = "Value") %>%
  dplyr::select("year", "partialRe", "Value") %>% 
  mutate(partialRe = sub("_.*", "", partialRe)) %>%
  filter(partialRe == "Re.juvv" | partialRe == "Re.aduu" | partialRe == "Re.bb" | 
           partialRe == "Re.resjuvv" | partialRe == "Re.resaduu" | partialRe == "Re.ress")

# Convert partialRe to factor with desired order
partialRe.year$partialRe <- factor(partialRe.year$partialRe, levels = c("Re.juvv", "Re.aduu", "Re.bb", "Re.resjuvv", "Re.resaduu", "Re.ress"))

partialRe.year <- partialRe.year %>%
  mutate(species = case_when(partialRe=="Re.juvv" ~ "Blackbird", 
                             partialRe=="Re.aduu" ~ "Blackbird", 
                             partialRe=="Re.bb" ~ "Blackbird", 
                             partialRe=="Re.resjuvv" ~ "Reservoir", 
                             partialRe=="Re.resaduu" ~ "Reservoir", 
                             partialRe=="Re.ress" ~ "Reservoir"))

plot.partialRe.route <-
  ggplot(data=partialRe.year) +
  facet_wrap(~species, scales = "free", nrow=2) +
  geom_violin(aes(x=factor(year), y=Value, fill=partialRe), position = position_dodge(width = 0.5)) + #, scale="width") +
  labs(x = "",
       y = "Partial Re",
       fill = "Transmission route") +
  scale_x_discrete(expand = expansion(mult = c(0.01, 0.01))) +
  scale_fill_manual(
    values = c("Re.juvv" = "#99FFCC", "Re.aduu" = "#66CC99", "Re.bb" = "#006600", "Re.resjuvv" = "#FFCCFF", "Re.resaduu" = "#FF99FF", "Re.ress" = "#CC6699"), 
    labels = c("Re.juvv" = "Juvenile blackbird - mosquito cycle",
               "Re.aduu" = "Adult blackbird - mosquito cycle",
               "Re.bb" = "Total blackbirds - mosquito cycle",
               "Re.resjuvv" = "Juvenile reservoir - mosquito cycle",
               "Re.resaduu" = "Adult reservoir - mosquito cycle", 
               "Re.ress" = "Total reservoir - mosquito cycle")) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.text=element_text(size=size), 
        legend.title=element_text(size=size), 
        # panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=size), 
        axis.title=element_text(size=size+2),
        strip.background = element_blank(),
        strip.text=element_text(size=size)) +
  guides(fill = guide_legend(ncol = 2))

plot.partialRe.route

ggsave(filename=paste0("../Output/Plots/dataFits/", runname, "/partialRe_routes_split.png", sep=""),
       width=13, height=7, bg="white")


# * R0 & Re together ---------------------------------------------------------
size <- 13
R.colors <- c("Re" = "coral",  "R0" = "aquamarine3")
R.colors <- c("Re" = "#666600",  "R0" = "#CCCC33")


plot.R0Re <-
  ggplot() +
  geom_line(data=Re.processed, aes(x=doy, y=mean.all, col="Re"), show.legend = FALSE) + 
  geom_line(data=R0.processed, aes(x=doy, y=mean.all, col="R0"), show.legend = FALSE) + 
  geom_ribbon(data=Re.processed, aes(x=doy , ymin=lb.all, ymax=ub.all, fill="Re"), alpha=0.5) +
  geom_ribbon(data=R0.processed, aes(x=doy , ymin=lb.all, ymax=ub.all, fill="R0"), alpha=0.5) +
  geom_text(data=meanRe, aes(x=110, y=15, label=paste("Mean Re: ", format(round(meanRe, digits=2), nsmall=2), sep="")), size=4) +
  geom_text(data=meanR0, aes(x=110, y=20, label=paste("Mean R0: ", format(round(meanR0, digits=2), nsmall=2), sep="")), size=4) +
  geom_hline(yintercept = 1, linetype="dashed") +
  
  facet_wrap(vars(year), nrow=7) +
  ylab("Reproduction number") + xlab('') +
  scale_x_continuous(breaks = c(92,153,214,275), 
                     labels = c("April","June","August","October")) +
  scale_color_manual(values=R.colors) +
  scale_fill_manual(values=R.colors) +
  
  theme_bw() +
  theme(legend.title=element_blank(), 
        legend.text = element_text(size=size),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=size), 
        axis.title=element_text(size=size+2),
        strip.background = element_blank(),
        strip.text=element_text(size=size))

plot.R0Re

ggsave(filename=paste0("../Output/Plots/dataFits/", runname, "/R0Re.png", sep=""), 
       width=6, height=11, bg="white")


# Live prevalence --------------------------------------------------------------

prev.processed <- process.prev(simOutput=simOutput, average="no")

prev.processed <- prev.processed %>%
  mutate(year = as.numeric(format(date,'%Y'))) %>%
  mutate(Week = week(date)) %>%
  mutate(doy = as.numeric(strftime(date, format = "%j"))) %>%
  mutate(month = as.numeric(format(date,'%m'))) %>%
  filter(month>=4 & month<=9)
  

prev.colors <- c("Live Juvenile" = "steelblue4", "Live Adult"= "steelblue1")

plot.liveprev <- 
  ggplot() +
  geom_line(data=prev.processed, aes(x=doy, y=mean.Juv, col="Live Juvenile"), show.legend = FALSE) +
  geom_line(data=prev.processed, aes(x=doy, y=mean.Adu, col="Live Adult"), show.legend = FALSE) +
  geom_ribbon(data=prev.processed, aes(x=doy , ymin=lb.Juv, ymax=ub.Juv, fill="Live Juvenile"), alpha=0.5) +
  geom_ribbon(data=prev.processed, aes(x=doy , ymin=lb.Adu, ymax=ub.Adu, fill="Live Adult"), alpha=0.5) +

  facet_wrap(vars(year), nrow=7) +
  ylim(0,0.1) +
  ylab("Live blackbird prevalence") + xlab("") +
  scale_x_continuous(breaks = c(92,153,214,275),
                     labels = c("April","June","August","October")) +
  
  scale_color_manual(values=prev.colors) +
  scale_fill_manual(values=prev.colors) +
  
  theme_bw() +
  theme(legend.title=element_blank(), 
        legend.text = element_text(size=size),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=size), 
        axis.title=element_text(size=size+2),
        strip.background = element_blank(),
        strip.text=element_text(size=size))


ggsave(filename=paste0("../Output/Plots/dataFits/", runname, "/liveprev_JuvAdu.png", sep="")
       , width=6.5, height=11, bg="white")


# Dead prevalence --------------------------------------------------------------

dead.processed <- process.dead(simOutput=simOutput, average="no")

dead.processed <- dead.processed %>%
  mutate(year = as.numeric(format(date,'%Y'))) %>%
  mutate(Week = week(date)) %>%
  mutate(doy = as.numeric(strftime(date, format = "%j"))) %>%
  mutate(month = as.numeric(format(date,'%m'))) %>%
  filter(month>=4 & month<=9)
  
prev.colors <- c("Dead Juvenile" = "firebrick" , "Dead Adult" = "darkorange")

plot.deadprev <-
  ggplot() +
  geom_line(data=dead.processed, aes(x=doy, y=mean.Juv, col="Dead Juvenile"), show.legend = FALSE) + 
  geom_line(data=dead.processed, aes(x=doy, y=mean.Adu, col="Dead Adult"), show.legend = FALSE) + 
  geom_ribbon(data=dead.processed, aes(x=doy , ymin=lb.Juv, ymax=ub.Juv, fill="Dead Juvenile"), alpha=0.5) +
  geom_ribbon(data=dead.processed, aes(x=doy , ymin=lb.Adu, ymax=ub.Adu, fill="Dead Adult"), alpha=0.5) +
  
  facet_wrap(vars(year), nrow=7) +
  ylim(0,1) +
  ylab("Dead blackbird prevalence") + xlab("") +
  scale_x_continuous(breaks = c(92,153,214,275),
                     labels = c("April","June","August","October")) +
  
  scale_color_manual(values=prev.colors) +
  scale_fill_manual(values=prev.colors) +
  
  theme_bw() +
  theme(legend.title=element_blank(), 
        legend.text = element_text(size=size),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=size), 
        axis.title=element_text(size=size+2),
        strip.background = element_blank(),
        strip.text=element_text(size=size))


ggsave(filename=paste0("../Output/Plots/dataFits/", runname, "/deadprev_JuvAdu.png", sep=""), 
       width=6.5, height=11, bg="white")



# Seroprevalence ----------------------------------------------------------

sero.processed <- process.sero(simOutput=simOutput, average="no")

sero.processed <- sero.processed %>%
  mutate(year = as.numeric(format(date,'%Y'))) %>%
  mutate(month = as.numeric(format(date,'%m'))) %>%
  mutate(Week = week(date)) %>%
  filter(month>=4 & month<=9) %>%
  mutate(doy = as.numeric(strftime(date, format = "%j")))

size <- 13
sero.colours <- c("Blackbird" = "mediumpurple", "Reservoir" = "goldenrod2" )

plot.sero <-
  ggplot() +
  geom_line(data=sero.processed, aes(x=doy, y=mean.bothbb, colour="Blackbird"),  show.legend = FALSE) + 
  geom_ribbon(data=sero.processed, aes(x=doy , ymin=lb.bothbb, ymax=ub.bothbb, fill="Blackbird"), alpha=0.5) +
  geom_line(data=sero.processed, aes(x=doy, y=mean.bothres, colour="Reservoir"),  show.legend = FALSE) +
  geom_ribbon(data=sero.processed, aes(x=doy , ymin=lb.bothres, ymax=ub.bothres, fill="Reservoir"), alpha=0.5) +
  
  facet_wrap(vars(year), nrow=7) +
  ylab("Seroprevalence") + xlab("") +
  scale_y_continuous(breaks = c(0,0.25,0.5, 0.75, 1)) + 
  # scale_x_continuous(breaks = c(1, 92,153,214,275), 
  #                    labels = c("January", "April","June","August","October")) +
  scale_x_continuous(breaks = c(1, 92, 183, 275), 
                     labels = c("January", "April","July", "October"),
                     expand = c(0, 0)) +
  coord_cartesian(ylim=c(0, 1)) +  
  scale_color_manual(values=sero.colours) +
  scale_fill_manual(values=sero.colours) +
  
  theme_bw() +
  theme(legend.title=element_blank(), 
        legend.text = element_text(size=size),
        legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=size), 
        axis.title=element_text(size=size+2),
        strip.background = element_blank(),
        strip.text=element_text(size=size))

ggsave(filename=paste0("../Output/Plots/dataFits/", runname, "/sero_both.png", sep=""),
       width=6, height=11, bg="white")



# Combine plots for main text panel ---------------------------------------

result.plotsB <- (plot.deadprev + plot.liveprev +  plot.R0Re) / plot.partialRe.route +
  plot_layout(heights = c(2, 1)) +   plot_annotation(tag_levels = 'A') 
result.plotsB

ggsave(filename=paste0("../Output/Plots/dataFits/", runname, "/panelB.png", sep=""), 
       plot=result.plotsB, 
       width=15, height=15, bg="white")
