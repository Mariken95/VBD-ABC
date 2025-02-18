
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

{library(SimInf)
library(tidyverse)
library(reshape2)     # to melt data into ggplot format
library(PerformanceAnalytics) # to make correlation plot
library(ggpubr)       # to arrange multiple ggplots
library(ggridges)     # for plots
library(philentropy)  # to calculate Kullback–Leibler distances
library(lubridate)    # faster as.Date transformations
library(patchwork)    # organise figure panel
library(RColorBrewer)
library(GGally)
library(HDInterval)}


# Functions combining chains ----------------------------------------------

# combine multiple posterior estimates and compare to input parameter

combine.chains <- function(abc.output) {
  particles <- read_csv(abc.output)
  # browser()
  
  particles <- particles %>%
    filter(gen == max(gen)) %>%
    dplyr::select(-c(pWeight, contains("dist") ))
  
  return(particles)
}

convergence.multichain <- function(abc.output) {
  particles <- read_csv(abc.output)
  # browser()
  
  particles <- particles %>%
    dplyr::select(-c(pWeight, contains("dist") ))
  
  return(particles)
  
}

extract.distance <- function(abc.output) {
  particles <- read_csv(abc.output)
  # browser()
  
  particles <- particles %>%
    filter(gen == max(gen)) %>%
    dplyr::select(contains("dist"))
  
  return(particles)
}

# Identifiability ------------------------------------------------------

runname <- "id_20240923"


# * Chosen values ---------------------------------------------------------

# set input parameters
parameters.A <- data.frame(matrix(nrow=7, ncol=3))
colnames(parameters.A) <- c("set", "Parameter", "SetValue")
parameters.A$Parameter <- c('vertical_transm', 'transmissionProbHM','scalingParameter', 
                          'disIndMortality', 'FOI.n', 'FOI.m', 'FOI.s') 
parameters.A$SetValue <- c(0.25, 0.7, 0.2, 0.2, 0.1, 0.25, 0.2) # set A
parameters.A$set <- "1"

parameters.B <- data.frame(matrix(nrow=7, ncol=3))
colnames(parameters.B) <- c("set", "Parameter", "SetValue")
parameters.B$Parameter <- c('vertical_transm', 'transmissionProbHM','scalingParameter', 
                            'disIndMortality', 'FOI.n', 'FOI.m', 'FOI.s')
parameters.B$SetValue <- c(0.25, 0.8, 0.15, 0.3, 0.08, 0.08, 0.08) # set B
parameters.B$set <- "2"

parameters.C <- data.frame(matrix(nrow=7, ncol=3))
colnames(parameters.C) <- c("set", "Parameter", "SetValue")
parameters.C$Parameter <- c('vertical_transm', 'transmissionProbHM','scalingParameter', 
                            'disIndMortality', 'FOI.n', 'FOI.m', 'FOI.s')
parameters.C$SetValue <- c(0.2, 0.9, 0.4, 0.2, 0.2, 0.2, 0.05) # set C
parameters.C$set <- "3"

parameters.D <- data.frame(matrix(nrow=7, ncol=3))
colnames(parameters.D) <- c("set", "Parameter", "SetValue")
parameters.D$Parameter <- c('vertical_transm', 'transmissionProbHM','scalingParameter', 
                            'disIndMortality', 'FOI.n', 'FOI.m', 'FOI.s')
parameters.D$SetValue <- c(0.2, 0.5, 0.15, 0.1, 0.04, 0.04, 0.04) # set D
parameters.D$set <- "4"

parameters.E <- data.frame(matrix(nrow=7, ncol=3))
colnames(parameters.E) <- c("set", "Parameter", "SetValue")
parameters.E$Parameter <- c('vertical_transm', 'transmissionProbHM','scalingParameter', 
                            'disIndMortality', 'FOI.n', 'FOI.m', 'FOI.s') 
parameters.E$SetValue <- c(0.1, 0.6, 0.1, 0.4, 0.01, 0.01, 0.1) # set E
parameters.E$set <- "5"

inputparameters <- rbind(parameters.A, parameters.B, parameters.C, parameters.D, parameters.E)
save(inputparameters, file="../Output/ABC/id_20240923/inputparameters.Rda")


# * Posteriors vs true ----------------------------------------------
load(paste0("../Output/ABC/", runname, "/inputparameters.Rda", sep=""))

# loop over chains and sets.
parameterSets <- c("1", "2", "3", "4", "5")
posterior.allchains <- data.frame()

for (set in parameterSets) {
  filename.temp <- paste("all_accepted_particles_", set, sep="")
  files.dir <- dir(paste0("../Output/ABC/", runname, "/", sep=""))
  files.temp <- files.dir[grep(filename.temp, files.dir, TRUE)]
  print(files.temp)
  
  for (chain in 1:length(files.temp)) {
    posterior <- combine.chains(abc.output=paste("../Output/ABC/", runname, "/", files.temp[chain], sep=""))
    print(head(posterior))
    posterior$chain <- chain
    posterior$set <- set
    posterior.allchains <- rbind(posterior.allchains, posterior)
    
  }
}

save(posterior.allchains, file=paste0("../Output/ABC/", runname, "/posterior.Rda", sep=""))

posteriors_long <- posterior.allchains %>%
    pivot_longer(cols=c(3:9), names_to = "Parameter", values_to = "Value")
posteriors_long <- left_join(posteriors_long, inputparameters)

ggplot(posteriors_long, aes(x=set, y=Value)) +
  facet_wrap(~ Parameter, ncol=4, scales = "free", labeller = labeller(Parameter = 
                                                                      c("disIndMortality" = "Disease-induced \nmortality rate blackbird", 
                                                                        "FOI.n" = "Historical FOI North",
                                                                        "FOI.m" = "Historical FOI Middle",
                                                                        "FOI.s" = "Historical FOI South",
                                                                        "scalingParameter" = "Abundance \nscaling parameter",
                                                                        "transmissionProbHM" = "Transmission probability \nhost-mosquito",
                                                                        "vertical_transm" = "Re-emergence \nrate"))) +
  geom_violin(scale="width", fill="lightblue", alpha=0.7, color=NA) +
  geom_point(aes(y=SetValue), shape=2, size=1.3) +
  stat_summary(fun="median", geom="point", col="darkblue", size=1.3, aes(shape="Estimated median")) +
  labs(x="Simulated parameter set", y="Value") +
  scale_shape_manual(name=NULL,
                     limits=c("Estimated median", "Simulated true value"),
                     values=c(19, 2)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), legend.position = c(0.87, 0.2))

ggsave(paste0("../Output/Plots/abcOutput/", runname, "/true-median.png"), width=8, height=4)

# * Coverage probability ----------------------------------------------------
# how many times is each true parameter included in the posterior
# compare true value to 95% HDI

load(paste0("../Output/ABC/", runname, "/inputparameters.Rda", sep=""))
load(paste0("../Output/ABC/", runname, "/posterior.Rda", sep=""))

get.coverage <- function(posterior, inputparameters, parameterset) {
  
  # browser()
  truevalue <- inputparameters %>% filter(set==parameterset)
  
  post.set <-  posterior %>%
    filter(set==set) %>% hdi(., credMass=0.95)
  post.set <- as.data.frame(post.set)
  post.set <- post.set %>%
    dplyr::select(-c(model, gen, chain, set))
  post.set <- t(post.set)
  post.set <- as.data.frame(post.set)
  post.set$true <- truevalue$SetValue
  post.set$contains <- ifelse(post.set$true<post.set$upper & post.set$true>post.set$lower, 1, 0)
  post.set$set <- parameterset
  
  return(post.set)
}

coverageA <- get.coverage(posterior=posterior.allchains, inputparameters=inputparameters, parameterset="1")
coverageB <- get.coverage(posterior=posterior.allchains, inputparameters=inputparameters, parameterset="2")
coverageC <- get.coverage(posterior=posterior.allchains, inputparameters=inputparameters, parameterset="3")
coverageD <- get.coverage(posterior=posterior.allchains, inputparameters=inputparameters, parameterset="4")
coverageE <- get.coverage(posterior=posterior.allchains, inputparameters=inputparameters, parameterset="5")

coverage <- rbind(coverageA, coverageB, coverageC, coverageD, coverageE)
coverage.group <- coverage %>% 
  group_by(set) %>%
  summarise(coverage.group = sum(contains)/7)

ggplot(data=coverage.group) +
  geom_col(aes(x=set, y=coverage.group)) +
  labs(x="Parameter set", y="Proportion of true values within 95% HDI") +
  ylim(0,1)

ggsave(filename=paste0("../Output/Plots/abcOutput/", runname, "/coverage.png", sep=""), bg="white")


# * Correlation true-estimate -----------------------------------------------
# pearson correlation coefficient between median and true value
load(paste0("../Output/ABC/", runname, "/inputparameters.Rda", sep=""))
load(paste0("../Output/ABC/", runname, "/posterior.Rda", sep=""))

get.correlation <- function(posterior, inputparameters, parameter) {
  
  browser()
  truevalue <- inputparameters %>% filter(Parameter==parameter)
  
  post.parameter <- posterior[c(parameter, 'set')] 
  
  post.parameter <- post.parameter %>%
    group_by(set) %>%
    mutate(row=row_number()) %>%
    pivot_wider(., names_from = set, values_from = paste(parameter))
  
  medians <- as.data.frame(apply(post.parameter[2:6], 2, median, na.rm=T))
  colnames(medians) <- c("median")
  medians$true <- truevalue$SetValue
  medians$correlation <- cor(x=medians$median, y=medians$true)
  medians$set <- c("1", "2", "3", "4", "5")
  medians$mean.true <- mean(medians$true)
  medians$mean.estim <- mean(medians$median)
  medians$sd.true <- sd(medians$true)
  medians$sd.estim <- sd(medians$median)
  medians$part1 <- (2*medians$sd.true*medians$sd.estim) / (medians$sd.true^2 + medians$sd.estim^2 + (medians$mean.true - medians$mean.estim)^2)
  medians$concordance <- medians$part1 * medians$correlation
  
  return(medians)
}

corr.vertical <- get.correlation(posterior=posterior.allchains, inputparameters=inputparameters, parameter='vertical_transm')
corr.tprob <- get.correlation(posterior=posterior.allchains, inputparameters=inputparameters, parameter='transmissionProbHM')
corr.scaling <- get.correlation(posterior=posterior.allchains, inputparameters=inputparameters, parameter='scalingParameter')
corr.mort <- get.correlation(posterior=posterior.allchains, inputparameters=inputparameters, parameter='disIndMortality')
corr.FOIn <- get.correlation(posterior=posterior.allchains, inputparameters=inputparameters, parameter='FOI.n')
corr.FOIm <- get.correlation(posterior=posterior.allchains, inputparameters=inputparameters, parameter='FOI.m')
corr.FOIs <- get.correlation(posterior=posterior.allchains, inputparameters=inputparameters, parameter='FOI.s')


corplot.vertical <- ggplot() +
  geom_point(data=corr.vertical, aes(x=median, y=true)) +
  geom_text(data=corr.vertical, aes(x=median, y=true, label=set), hjust=0, vjust=0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  annotate("text", x=0.09, y=0.25, label=paste("Concordance: ", round(corr.vertical$concordance[1], 2), sep=""), size=5) +
  labs(x="Estimated median", y="True value", title = "Re-emergence rate") +
  xlim(0, 0.3) + ylim(0, 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.text=element_text(size=16), axis.title=element_text(size=16),
        plot.title = element_text(size = 16, face = "bold"))

corplot.transm <- ggplot() +
  geom_point(data=corr.tprob, aes(x=median, y=true)) +
  geom_text(data=corr.tprob, aes(x=median, y=true, label=set), hjust=0, vjust=0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  annotate("text", x=0.23, y=0.75, label=paste("Concordance: ", round(corr.tprob$concordance[1], 2), sep=""), size=5) +
  labs(x="Estimated median", y="True value", title = "Transmission probability \nhost-mosquito") +
  xlim(0, 0.9) + ylim(0, 0.9) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text=element_text(size=16), axis.title=element_text(size=16),
        plot.title = element_text(size = 16, face = "bold"))

corplot.scaling <- ggplot() +
  geom_point(data=corr.scaling, aes(x=median, y=true)) +
  geom_text(data=corr.scaling, aes(x=median, y=true, label=set), hjust=0, vjust=0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  annotate("text", x=0.27, y=0.82, label=paste("Concordance: ", round(corr.scaling$concordance[1], 2), sep=""), size=5) +
  labs(x="Estimated median", y="True value", title= "Abundance scaling parameter") +
  xlim(0, 1) + ylim(0, 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text=element_text(size=16), axis.title=element_text(size=16),
        plot.title = element_text(size = 16, face = "bold"))

corplot.mort <- ggplot() +
  geom_point(data=corr.mort, aes(x=median, y=true)) +
  geom_text(data=corr.mort, aes(x=median, y=true, label=set), hjust=0, vjust=0) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  annotate("text", x=0.15, y=0.5, label=paste("Concordance: ", round(corr.mort$concordance[1], 2), sep=""), size=5) +
  labs(x="Estimated median", y="True value", title="Disease-induced \nmortality rate blackbird") +
  xlim(0, 0.6) + ylim(0, 0.6) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text=element_text(size=16), axis.title=element_text(size=16),
        plot.title = element_text(size = 16, face = "bold"))

corplot.FOI <- ggplot() +
  geom_point(data=corr.FOIn, aes(x=median, y=true, colour="Historical FOI north")) +
  annotate("text", x=0.05, y=0.26, label=paste("North: ", round(corr.FOIn$concordance[1], 2), sep=""), size=5) +
  
  geom_point(data=corr.FOIm, aes(x=median, y=true, colour="Historical FOI middle")) +
  annotate("text", x=0.05, y=0.23, label=paste("Middle: ", round(corr.FOIm$concordance[1], 2), sep=""), size=5) +
  
  geom_point(data=corr.FOIs, aes(x=median, y=true, colour="Historical FOI south")) +
  annotate("text", x=0.05, y=0.20, label=paste("South: ", round(corr.FOIs$concordance[1], 2), sep=""), size=5) +
  
  labs(x="Estimated median", y="True value", colour="Parameter", title="Historical FOI") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  
  xlim(0, 0.3) + ylim(0, 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text=element_text(size=16), axis.title=element_text(size=16),
        plot.title = element_text(size = 16, face = "bold"))

corplot.FOI.bw <- 
  ggplot() +
  geom_point(data=corr.FOIn, aes(x=median, y=true, shape="Historical FOI north")) +
  annotate("text", x=0.05, y=0.26, label=paste("North: ", round(corr.FOIn$concordance[1], 2), sep=""), size=5) +
  
  geom_point(data=corr.FOIm, aes(x=median, y=true, shape="Historical FOI middle")) +
  annotate("text", x=0.05, y=0.23, label=paste("Middle: ", round(corr.FOIm$concordance[1], 2), sep=""), size=5) +
  
  geom_point(data=corr.FOIs, aes(x=median, y=true, shape="Historical FOI south")) +
  annotate("text", x=0.05, y=0.20, label=paste("South: ", round(corr.FOIs$concordance[1], 2), sep=""), size=5) +
  
  labs(x="Estimated median", y="True value", shape="Parameter", title="Historical FOI") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  
  xlim(0, 0.3) + ylim(0, 0.3) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text=element_text(size=16), axis.title=element_text(size=16), 
        legend.text=element_text(size=14), legend.title=element_text(size=14), 
        plot.title = element_text(size = 16, face = "bold"))


(corplot.vertical / corplot.transm) | (corplot.scaling / corplot.mort) | (corplot.FOI / plot_spacer() )
ggsave(filename=paste0("../Output/Plots/abcOutput/", runname, "/correlation_truemedian.png", sep=""),bg="white", width=14, height=8)

(corplot.vertical / corplot.transm) | (corplot.scaling / corplot.mort) | (corplot.detection / corplot.FOI.bw)
ggsave(filename="../Output/Plots/abcOutput/runs_20231031/correlation_truemedian_bw.png", bg="white", width=14.5, height=8)


# Estimation from data ----------------------------------------------------

runname <- "runs_20240805"
version <- "D"

# * Posteriors ------------------------------------------------------------

# loop over chains and sets
posterior.allchains <- data.frame()

filename.temp <- paste0("accepted_particles_", version, sep="")
files.dir <- dir(paste0("../Output/ABC/", runname, "/", sep=""))
files.temp <- files.dir[grep(filename.temp, files.dir, TRUE)]
for (chain in 1:length(files.temp)) {
  posterior <- combine.chains(abc.output=paste("../Output/ABC/", runname, "/", files.temp[chain], sep=""))
  print(head(posterior))
  posterior$chain <- chain
  posterior.allchains <- rbind(posterior.allchains, posterior)
}

save(posterior.allchains, file=paste0("../Output/ABC/", runname, "/posterior", version, ".Rda", sep=""))

posteriors_long <- posterior.allchains %>%
  pivot_longer(cols=c(3:9), names_to = "Parameter", values_to = "Value")

# ** values of posterior with CI ------------------------------------------

## overwintering
median(posterior.allchains$vertical_transm)
hdi(posterior.allchains$vertical_transm, credMass=0.95)

## transmission prob
median(posterior.allchains$transmissionProbHM)
hdi(posterior.allchains$transmissionProbHM, credMass=0.95)

## scaling parameter
median(posterior.allchains$scalingParameter)
hdi(posterior.allchains$scalingParameter, credMass=0.95)

## disease-induced mortality / recovery rate
median(posterior.allchains$disIndMortality)
hdi(posterior.allchains$disIndMortality, credMass=0.95)

median(posterior.allchains$disIndMortality) / (median(posterior.allchains$disIndMortality) + 0.25)
dis.mort <- hdi(posterior.allchains$disIndMortality, credMass=0.95)
dis.mort[1]/ (dis.mort[1]+0.25)
dis.mort[2]/ (dis.mort[2]+0.25)

## historical FOI
median(posterior.allchains$FOI)
hdi(posterior.allchains$FOI, credMass=0.95)

## lifespan
median(posterior.allchains$lifespanR)/365.25
hdi(posterior.allchains$lifespanR, credMass=0.95)/365
# lifespan at birth. juv mortality - 0.0059. adult mortality is 1/lifespan
# lifespan at birth = lifespan juv + P(survive as juv) * lifespan adult
# lifespan juv = (1 - (1-0.0059)^345) / 0.0059
# P(survive as juv) = (1-0.0059)^345
# lifespan adult = median(posterior.allchains$lifespanR)
( (1 - (1-0.0059)^345) / 0.0059 + (1-0.0059)^345 * median(posterior.allchains$lifespanR) ) / 365.25 # 1.3 years

# for blackbird
( (1 - (1-0.0059)^345) / 0.0059 + (1-0.0059)^345 * (1/0.0011) ) / 365.25 # 0.72 years

## biting pref
median(posterior.allchains$bitingPref)
hdi(posterior.allchains$bitingPref, credMass=0.95)


# * Convergence -----------------------------------------------------------

filename.temp <- "accepted_particles_D"
files.dir <- dir(paste0("../Output/ABC/", runname, "/", sep=""))
files.temp <- files.dir[grep(filename.temp, files.dir, TRUE)]

convergence.allchains <- data.frame()

for (chain in 1:length(files.temp)) {
  convergence <- convergence.multichain(abc.output=paste("../Output/ABC/", runname, "/", files.temp[chain], sep=""))
  print(head(convergence))
  convergence$chain <- chain
  convergence.allchains <- rbind(convergence.allchains, convergence)
}


convergence_long <- convergence.allchains %>%   
  pivot_longer(cols=c(3:9), names_to = "Parameter", values_to = "Value")

ggplot(convergence_long, aes(x=Value,y=gen, group=gen))+
  facet_wrap(~Parameter, scales = "free", ncol=4, labeller = labeller(Parameter = 
                                                                        c("bitingPref" = "Relative biting \nfrequency parameter ", 
                                                                          "disIndMortality" = "Disease-induced \nmortality rate blackbird", 
                                                                          "FOI" = "Historical FOI",
                                                                          "lifespanR" = "Lifespan \nreservoir population",
                                                                          "scalingParameter" = "Abundance \nscaling parameter",
                                                                          "transmissionProbHM" = "Transmission probability \nhost-mosquito",
                                                                          "vertical_transm" = "Re-emergence \nrate"))) +
  geom_density_ridges(aes(), scale = 1.6, alpha = .6, color = "#aaaaaa00") + theme_ridges()+
  geom_text(aes(x=-0.2, y=gen, label=gen)) +
  scale_y_discrete(expand = c(0.02, 0)) +
  labs(x="Value", y="Generation")

ggsave(paste0("../Output/Plots/abcOutput/", runname, "/convergence_D.png", sep=""), width=12, height=8, bg = "white")


# * Correlation -------------------------------------------------------------

colnames(posterior.allchains) <- c("gen", "model", "Re-emergence \nrate", "Transmission probability \nhost-mosquito",
                                   "Abundance \nscaling parameter", "Disease-induced \nmortality rate blackbird",
                                   "Historical FOI", "Lifespan \nreservoir population", "Relative biting \nfrequency parameter", 
                                   "chain")
chart.Correlation(posterior.allchains[,3:9], histogram = TRUE, method = "pearson")
cor(posterior.allchains[,3:7])



# * Model comparison ------------------------------------------------------
runname <- "runs_20240805"
# A
distances.allchains <- data.frame()

filename.temp <- "accepted_particles_A"
files.dir <- dir(paste0("../Output/ABC/", runname, "/", sep=""))
files.temp <- files.dir[grep(filename.temp, files.dir, TRUE)]
for (chain in 1:length(files.temp)) {
  distance <- extract.distance(abc.output=paste("../Output/ABC/", runname, "/", files.temp[chain], sep=""))
  print(head(distance))
  distance$chain <- chain
  distances.allchains <- rbind(distances.allchains, distance)
}

distance.A <- distances.allchains 

# B. 
distances.allchains <- data.frame()

filename.temp <- "accepted_particles_B"
files.dir <- dir(paste0("../Output/ABC/", runname, "/", sep=""))
files.temp <- files.dir[grep(filename.temp, files.dir, TRUE)]
for (chain in 1:length(files.temp)) {
  distance <- extract.distance(abc.output=paste("../Output/ABC/", runname, "/", files.temp[chain], sep=""))
  print(head(distance))
  distance$chain <- chain
  distances.allchains <- rbind(distances.allchains, distance)
}

distance.B <- distances.allchains 

# C.
distances.allchains <- data.frame()

filename.temp <- "accepted_particles_C"
files.dir <- dir(paste0("../Output/ABC/", runname, "/", sep=""))
files.temp <- files.dir[grep(filename.temp, files.dir, TRUE)]
for (chain in 1:length(files.temp)) {
  distance <- extract.distance(abc.output=paste("../Output/ABC/", runname, "/", files.temp[chain], sep=""))
  print(head(distance))
  distance$chain <- chain
  distances.allchains <- rbind(distances.allchains, distance)
}

distance.C <- distances.allchains 

# D 
distances.allchains <- data.frame()

filename.temp <- "accepted_particles_D"
files.dir <- dir(paste0("../Output/ABC/", runname, "/", sep=""))
files.temp <- files.dir[grep(filename.temp, files.dir, TRUE)]
for (chain in 1:length(files.temp)) {
  distance <- extract.distance(abc.output=paste("../Output/ABC/", runname, "/", files.temp[chain], sep=""))
  print(head(distance))
  distance$chain <- chain
  distances.allchains <- rbind(distances.allchains, distance)
}

distance.D <- distances.allchains 

# S1
distances.allchains <- data.frame()

filename.temp <- "accepted_particles_S1"
files.dir <- dir(paste0("../Output/ABC/", runname, "/", sep=""))
files.temp <- files.dir[grep(filename.temp, files.dir, TRUE)]
for (chain in 1:length(files.temp)) {
  distance <- extract.distance(abc.output=paste("../Output/ABC/", runname, "/", files.temp[chain], sep=""))
  print(head(distance))
  distance$chain <- chain
  distances.allchains <- rbind(distances.allchains, distance)
}

distance.S1 <- distances.allchains 

# S2
distances.allchains <- data.frame()

filename.temp <- "accepted_particles_S2"
files.dir <- dir(paste0("../Output/ABC/", runname, "/", sep=""))
files.temp <- files.dir[grep(filename.temp, files.dir, TRUE)]
for (chain in 1:length(files.temp)) {
  distance <- extract.distance(abc.output=paste("../Output/ABC/", runname, "/", files.temp[chain], sep=""))
  print(head(distance))
  distance$chain <- chain
  distances.allchains <- rbind(distances.allchains, distance)
}

distance.S2 <- distances.allchains 

# S3
distances.allchains <- data.frame()

filename.temp <- "accepted_particles_S3"
files.dir <- dir(paste0("../Output/ABC/", runname, "/", sep=""))
files.temp <- files.dir[grep(filename.temp, files.dir, TRUE)]
for (chain in 1:length(files.temp)) {
  distance <- extract.distance(abc.output=paste("../Output/ABC/", runname, "/", files.temp[chain], sep=""))
  print(head(distance))
  distance$chain <- chain
  distances.allchains <- rbind(distances.allchains, distance)
}

distance.S3 <- distances.allchains 

## comparison
distances.comparison <- rbind(colMeans(distance.A),
                              colMeans(distance.B),
                              colMeans(distance.C),
                              colMeans(distance.D),
                              colMeans(distance.S1),
                              colMeans(distance.S2),
                              colMeans(distance.S3))

distances.comparison.norm <- as.data.frame(apply(distances.comparison[,1:5], 2, function(x){x/max(x)}))
colnames(distances.comparison.norm) <- c("Untested \n dead bird", "Prevalence \n live birds", 
                                         "Seroprevalence \n live birds", "Prevalence \n dead birds", "Population size")

distances.comparison.norm$scenario <- rep(c("A: Blackbird only", 
                                            "B: Reservoir host: increased dispersal", 
                                            "C: Reservoir host: no infection mortality & \n    estimated lifespan",
                                            "D: Reservoir host: increased dispersal & \n    no infection mortality & estimated lifespan",
                                            "S1: Introductions country-wide",
                                            "S2: Uniform spatial distribution reservoir host",
                                            "S3: Alternative movement distances"),
                                          each=1)

distances.comparison_long <- distances.comparison.norm %>%
  pivot_longer(cols=c(1:5), names_to = "Distance", values_to = "Value")
distances.comparison_long$scenario <- factor(distances.comparison_long$scenario, 
                                             levels = c("A: Blackbird only", 
                                                        "B: Reservoir host: increased dispersal", 
                                                        "C: Reservoir host: no infection mortality & \n    estimated lifespan",
                                                        "D: Reservoir host: increased dispersal & \n    no infection mortality & estimated lifespan",
                                                        "S1: Introductions country-wide",
                                                        "S2: Uniform spatial distribution reservoir host",
                                                        "S3: Alternative movement distances"))

## main scenarios
distances.comparison_long.main <- distances.comparison_long[1:20,]

ggplot(distances.comparison_long.main, aes(x=Distance, y=Value, fill=scenario)) +
  geom_bar(stat="identity", position = position_dodge2(reverse = TRUE)) +
  ylab("Normalised distance to data, \n relative to max distance") + xlab("") +
  scale_fill_manual(values = brewer.pal(n = 8, name = 'Greens')[c(3, 4, 6, 7)]) +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text=element_text(size=12), legend.text=element_text(size=12), legend.title = element_blank())

ggsave("../Output/Plots/abcOutput/comparisonsSep/distance_comparison_main2.png", width=11, height=5)

## sensitivity scenarios
distances.comparison_long.S <- distances.comparison_long[16:35,]

ggplot(distances.comparison_long.S, aes(x=Distance, y=Value, fill=scenario)) +
  geom_bar(stat="identity", position = position_dodge2(reverse = TRUE)) +
  ylab("Normalised distance to data, \n relative to max distance") + xlab("") +
  scale_fill_manual(values = c(brewer.pal(n = 8, name = 'Greens')[c(7)], brewer.pal(n = 8, name = 'Blues')[c(3, 5, 7)])) +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text=element_text(size=12), legend.text=element_text(size=12), legend.title = element_blank())

ggsave("../Output/Plots/abcOutput/comparisonsSep/distance_comparison_S_2.png", width=11, height=5)

rowMeans(distances.comparison.norm[,1:5])

# * Acceptance rate ----------------------------------------------
threshold.allchains <- data.frame()

filename.temp <- "thresholds_D"
files.dir <- dir(paste0("../Output/ABC/", runname, "/", sep=""))
files.temp <- files.dir[grep(filename.temp, files.dir, TRUE)]

for (chain in 1:length(files.temp)) {
  threshold <- read_csv(paste("../Output/ABC/", runname, "/", files.temp[chain], sep=""))
  print(head(threshold))
  threshold$chain <- chain
  threshold.allchains <- rbind(threshold.allchains, threshold)
}

threshold.allchains$chain <- as.factor(threshold.allchains$chain)

ggplot(threshold.allchains, aes(x=as.factor(gen), y=acceptance_rate)) +
  geom_point(aes(y=acceptance_rate, colour=chain), size=2 ) +
  geom_hline(yintercept = 0.05, linetype="dashed") +
  xlab("Generation") + ylab("Acceptance rate") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(paste0("../Output/Plots/abcOutput/", runname, "/acceptance_D.png", sep=""), width=5, height=4)


# * Distances for each generation  ----------------------------------------------
thresholds_long <- threshold.allchains %>%
  dplyr::select(-c(acceptance_rate)) %>%
  pivot_longer(cols=c(2:6), names_to = "Pattern", values_to = "Threshold")

ggplot(thresholds_long, aes(x=as.factor(gen), y=Threshold)) +
  facet_wrap(~Pattern, scales = "free", ncol=5, labeller = labeller(Pattern = 
                                                                      c("dist1" = "Untested \n dead blackbird", 
                                                                        "dist2" = "Prevalence \n live blackbirds", 
                                                                        "dist3" = "Seroprevalence \n live blackbirds",
                                                                        "dist4" = "Prevalence \n dead blackbirds",
                                                                        "dist5" = "Relative \n population size"))) +
  geom_point(aes(y=Threshold), shape=2, size=2) +
  xlab("Generation") + ylab("Acceptance threshold value") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

ggsave(paste0("../Output/Plots/abcOutput/", runname, "/thresholdsByGen_D.png", sep=""), width=8, height=3)
