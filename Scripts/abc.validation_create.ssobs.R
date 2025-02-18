# This script creates the observed data patterns from specified parameter sets, for identifiability runs

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(SimInf)
library(tidyverse)
library(lubridate)    
library(data.table)
library(patchwork)

# load info needed to calculate summary statistics
source("abc_functions_summaryStatistics.R")                     # loads functions to calculate summary statistics for ABC

load("../Data/abc_summaryStatistics/PCRprev_likelihood.Rda")           
load("../Data/abc_summaryStatistics/sero_likelihood.Rda")      
load("../Data/abc_summaryStatistics/dead.prev_likelihood.Rda")          

# load simulated trajectories & calculate summary statistics
filename.temp <- "traject"
files.temp <- dir("../Data/abc_identifiability/")
files.temp <- files.temp[grep(filename.temp, files.temp, TRUE)]

ss_obs <- list()

for (set in 1:length(files.temp)) {
  print(set)
  load(paste0("../Data/abc_identifiability/", files.temp[set]))
  traject.withindex <- prep.modeloutput(modeloutput=traject)
  ss_sim1 <- calculate.deadbird_year.pattern(traject.withindex=traject.withindex, output.type="statistic") 
  ss_sim2 <- calculate.prevalence.likelihood(traject.withindex=traject.withindex, PCRdataset=PCRprev_likelihood, output.type="statistic")
  ss_sim3 <- calculate.seroprev.likelihood(traject.withindex=traject.withindex, serodataset=sero_likelihood, output.type="statistic")
  ss_sim4 <- calculate.deadprevalence.likelihood(traject.withindex=traject.withindex, deaddataset=dead.prev_likelihood, output.type="statistic")
  ss_sim5 <- calculate.birdpop(traject.withindex=traject.withindex, output.type="statistic")

  ss_sim <- list(ss_sim1, ss_sim2, ss_sim3, ss_sim4, ss_sim5)
  
  ss_obs[[set]] <- ss_sim
}


save(ss_obs, file="../Data/abc_identifiability/summStat.Rda")


# Plots of summary statistics ---------------------------------------------

# * set A -----------------------------------------------------------------
load("../Data/abc_identifiability/traject_A.Rda")   

A.plot <- prep.modeloutput(modeloutput=traject)
A.deadbird <-  calculate.deadbird_year.pattern(traject.withindex=A.plot, output.type="plot")
A.liveprev <- calculate.prevalence.likelihood(traject.withindex=A.plot, PCRdataset=PCRprev_likelihood, output.type="plot")
A.seroprev <- calculate.seroprev.likelihood(traject.withindex=A.plot, serodataset=sero_likelihood, output.type="plot")
A.deadprev <- calculate.deadprevalence.likelihood(traject.withindex=A.plot, deaddataset=dead.prev_likelihood, output.type="plot")
A.pop <- calculate.birdpop(traject.withindex=A.plot, output.type="plot")


pattern_plots <- (A.deadbird + A.pop) / (A.liveprev + A.seroprev + A.deadprev) + 
  plot_annotation(tag_levels = 'A')
pattern_plots

ggsave(filename="../Data/abc_identifiability/summaryStatistics_A.png", plot=pattern_plots, 
       width=13, height=10, bg="white")


# * set B -------------------------------------------------------------------
load("../Data/abc_identifiability/traject_B.Rda")                  
B.plot <- prep.modeloutput(modeloutput=traject)
B.deadbird <-  calculate.deadbird_year.pattern(traject.withindex=B.plot, output.type="plot")
B.liveprev <- calculate.prevalence.likelihood(traject.withindex=B.plot, PCRdataset=PCRprev_likelihood, output.type="plot")
B.seroprev <- calculate.seroprev.likelihood(traject.withindex=B.plot, serodataset=sero_likelihood, output.type="plot")
B.deadprev <- calculate.deadprevalence.likelihood(traject.withindex=B.plot, deaddataset=dead.prev_likelihood, output.type="plot")
B.pop <- calculate.birdpop(traject.withindex=B.plot, output.type="plot")

pattern_plots <- (B.deadbird + B.pop) / (B.liveprev + B.seroprev + B.deadprev) + 
  plot_annotation(tag_levels = 'A')
pattern_plots

ggsave(filename="../Data/abc_identifiability/summaryStatistics_B.png", plot=pattern_plots, 
       width=13, height=10, bg="white")

# * set C -------------------------------------------------------------------
load("../Data/abc_identifiability/traject_C.Rda")                  
C.plot <- prep.modeloutput(modeloutput=traject)
C.deadbird <-  calculate.deadbird_year.pattern(traject.withindex=C.plot, output.type="plot")
C.liveprev <- calculate.prevalence.likelihood(traject.withindex=C.plot, PCRdataset=PCRprev_likelihood, output.type="plot")
C.seroprev <- calculate.seroprev.likelihood(traject.withindex=C.plot, serodataset=sero_likelihood, output.type="plot")
C.deadprev <- calculate.deadprevalence.likelihood(traject.withindex=C.plot, deaddataset=dead.prev_likelihood, output.type="plot")
C.pop <- calculate.birdpop(traject.withindex=C.plot, output.type="plot")

pattern_plots <- (C.deadbird + C.pop) / (C.liveprev + C.seroprev + C.deadprev) + 
  plot_annotation(tag_levels = 'A')
pattern_plots

ggsave(filename="../Data/abc_identifiability/summaryStatistics_C.png", plot=pattern_plots, 
       width=13, height=10, bg="white")

# * set D -------------------------------------------------------------------
load("../Data/abc_identifiability/traject_D.Rda")                  
D.plot <- prep.modeloutput(modeloutput=traject)
D.deadbird <-  calculate.deadbird_year.pattern(traject.withindex=D.plot, output.type="plot")
D.liveprev <- calculate.prevalence.likelihood(traject.withindex=D.plot, PCRdataset=PCRprev_likelihood, output.type="plot")
D.seroprev <- calculate.seroprev.likelihood(traject.withindex=D.plot, serodataset=sero_likelihood, output.type="plot")
D.deadprev <- calculate.deadprevalence.likelihood(traject.withindex=D.plot, deaddataset=dead.prev_likelihood, output.type="plot")
D.pop <- calculate.birdpop(traject.withindex=D.plot, output.type="plot")

pattern_plots <- (D.deadbird + D.pop) / (D.liveprev + D.seroprev + D.deadprev) + 
  plot_annotation(tag_levels = 'A')
pattern_plots

ggsave(filename="../Data/abc_identifiability/summaryStatistics_D.png", plot=pattern_plots, 
       width=13, height=10, bg="white")

# * set E -------------------------------------------------------------------
load("../Data/abc_identifiability/traject_E.Rda")                  
E.plot <- prep.modeloutput(modeloutput=traject)
E.deadbird <-  calculate.deadbird_year.pattern(traject.withindex=E.plot, output.type="plot")
E.liveprev <- calculate.prevalence.likelihood(traject.withindex=E.plot, PCRdataset=PCRprev_likelihood, output.type="plot")
E.seroprev <- calculate.seroprev.likelihood(traject.withindex=E.plot, serodataset=sero_likelihood, output.type="plot")
E.deadprev <- calculate.deadprevalence.likelihood(traject.withindex=E.plot, deaddataset=dead.prev_likelihood, output.type="plot")
E.pop <- calculate.birdpop(traject.withindex=E.plot, output.type="plot")

pattern_plots <- (E.deadbird + E.pop) / (E.liveprev + E.seroprev + E.deadprev) + 
  plot_annotation(tag_levels = 'A')
pattern_plots

ggsave(filename="../Data/abc_identifiability/summaryStatistics_E.png", plot=pattern_plots, 
       width=13, height=10, bg="white")

