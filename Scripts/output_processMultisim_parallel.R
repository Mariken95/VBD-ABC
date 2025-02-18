# PROCESS MULTIPLE SIMULATION OUTPUT FOR DATA FIT & R0

# abc_diagnostics saves posterior samples
# main_run_multisim.R runs repeated simulations with posterior samples
# here: calculate summ stat per simulation & bring together for plots 


# HPC ------------------------------------------------------------------
rm(list=ls())

# test 
print('is something happening?')

# load libraries -

.libPaths( c("/home/mdwit/R/library" , .libPaths() ) )

library(tidyverse,lib="/home/mdwit/R/library")
library(data.table,lib="/home/mdwit/R/library")
library(zoo,lib="/home/mdwit/R/library")
library(optparse,lib="/home/mdwit/R/library")

# load files to calculate summary statistics
load("../Data/abc_summaryStatistics/PCRprev_likelihood.Rda")           
load("../Data/abc_summaryStatistics/sero_likelihood.Rda")      
load("../Data/abc_summaryStatistics/dead.prev_likelihood.Rda")          

source("functions_processOutput.R") # for all non-summary statistics calculations
source("abc_functions_summaryStatistics.R")


# set-up ----------------------------------------------------------

# Define the command-line arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL,
              help = "input file path", metavar = "character"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "output file path", metavar = "character")
)

# Parse the command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if input file was provided
if (is.null(opt$input)) {
  stop("Input file must be provided")
}

# Check if output file was provided
if (is.null(opt$output)) {
  stop("Output file must be provided")
}

print(opt$input)
print(opt$output)


# Create output folder
folder_path <- opt$output
# Check if the folder path exists
if (!dir.exists(folder_path)) {
  # Folder does not exist, create the folder
  if (dir.create(folder_path, recursive = TRUE)) {
   cat("Folder created successfully.\n")
  } else {
   cat("Error: Unable to create folder.\n")
  }
} else {
  cat("Folder already exists.\n")
}

# Read the input file (called 'listsofresults')
load(opt$input)


# Summary statistics from simulations ------------------------------

ss_sim_posterior <- list()

for (run in 1:length(listofresults)) {

  print(run)

  traject.withindex <- prep.modeloutput(modeloutput=listofresults[[run]][[1]])
  ss_sim1 <- calculate.deadbird_year.pattern(traject.withindex=traject.withindex, output.type="table")
  ss_sim2 <- calculate.prevalence.likelihood(traject.withindex=traject.withindex, PCRdataset=PCRprev_likelihood, output.type="table")
  ss_sim3 <- calculate.seroprev.likelihood(traject.withindex=traject.withindex, serodataset=sero_likelihood, output.type="table")
  ss_sim4 <- calculate.deadprevalence.likelihood(traject.withindex=traject.withindex, deaddataset=dead.prev_likelihood, output.type="table")
  ss_sim5 <- calculate.birdpop(traject.withindex=traject.withindex, output.type="table")

  ss_sim <- list(ss_sim1, ss_sim2, ss_sim3, ss_sim4, ss_sim5)
  ss_sim_posterior[[run]] <- ss_sim
}

saveRDS(ss_sim_posterior, file=sprintf(paste(opt$output, '/summ_stat_posterior.rds', sep='')))


# R and prevalence from simulations ---------------------------------------------------------------

# * Time: Calculate R0, Re, prevalence --------------------------------------------

simOutput <- list()

for (run in 1:length(listofresults)) { 

  print(run)
  R0_sim <- list()
  prop.S <- list()
  Re_sim <- list()
  prev <- list()
  sero <- list()
  dead <- list()
  biting <- list()

  modeloutput <- listofresults[[run]][['modeloutput']]

  modeloutput <- modeloutput %>%
    mutate(date = as.Date(time, origin="2015-12-31")) %>%
    mutate(month = as.numeric(format(date,'%m'))) %>%
    mutate(year = as.numeric(format(date,'%Y'))) #%>%

  # aggregate across all locations
  setDT(modeloutput)
  modeloutput <- rollup(modeloutput, j=lapply(.SD, mean), by=c("date"))
  # this function adds a final row with the colMeans -> needs to be removed
  modeloutput <- modeloutput[1:(nrow(modeloutput)-1),]

  traject.withindex <- prep.modeloutput(modeloutput=modeloutput)

  traject.withindex <- calculate.newdeaths(output=traject.withindex)

  for(row in 1:nrow(traject.withindex)) {
    browser()
    R0_sim_row <- calculate.R0(output=traject.withindex[row,], posterior=listofresults[[run]][[2]],
                               transmissionProbMH=0.88, IIR=0.67, mortality.adu=0.0011, mortality.juv=0.0059, recovery=0.25,
                               version="multi-host")
    R0_sim <- rbind(R0_sim, R0_sim_row)
    
    prop.S_row <- calculate.propS(output=traject.withindex[row,])
    prop.S <- rbind(prop.S, prop.S_row)
    
    Re_sim_row <- calculate.Re(R0_sim = R0_sim[row,], prop.S = prop.S[row,])
    Re_sim <- rbind(Re_sim, Re_sim_row)
    
    prev_row <- calculate.prev(output=traject.withindex[row,])
    prev <- rbind(prev, prev_row)

    sero_row <- calculate.sero(output=traject.withindex[row,])
    sero <- rbind(sero, sero_row)

    dead_row <- calculate.dead(output=traject.withindex[row,])
    dead <- rbind(dead, dead_row)

    biting_row <- calculate.biting.bb(output=traject.withindex[row,], posterior=listofresults[[run]][[2]])
    biting <- rbind(biting, biting_row)

  }
  traject.withindex$R0.mosq.juv <- abs(as.numeric(R0_sim[,1]))
  traject.withindex$R0.mosq.adu <- abs(as.numeric(R0_sim[,2]))
  traject.withindex$R0.mosq.resjuv <- abs(as.numeric(R0_sim[,3]))
  traject.withindex$R0.mosq.resadu <- abs(as.numeric(R0_sim[,4]))
  traject.withindex$R0.juv.mosq <- abs(as.numeric(R0_sim[,5]))
  traject.withindex$R0.adu.mosq <- abs(as.numeric(R0_sim[,6]))
  traject.withindex$R0.resjuv.mosq <- abs(as.numeric(R0_sim[,7]))
  traject.withindex$R0.resadu.mosq <- abs(as.numeric(R0_sim[,8]))
  traject.withindex$R0.juv <- abs(as.numeric(R0_sim[,9]))
  traject.withindex$R0.adu <- abs(as.numeric(R0_sim[,10]))
  traject.withindex$R0.bb <- abs(as.numeric(R0_sim[,11]))
  traject.withindex$R0.resjuv <- abs(as.numeric(R0_sim[,12]))
  traject.withindex$R0.resadu <- abs(as.numeric(R0_sim[,13]))
  traject.withindex$R0.res <- abs(as.numeric(R0_sim[,14]))
  traject.withindex$R0 <- abs(as.numeric(R0_sim[,15]))
  traject.withindex$R0_7d <- rollmeanr(traject.withindex$R0, 7, fill=NA)

  traject.withindex$prop.sus.juv <- as.numeric(prop.S[,1])
  traject.withindex$prop.sus.adu <- as.numeric(prop.S[,2])
  traject.withindex$prop.sus.resjuv <- as.numeric(prop.S[,3])
  traject.withindex$prop.sus.resadu <- as.numeric(prop.S[,4])
  traject.withindex$prop.sus.mosq <- as.numeric(prop.S[,5])

  traject.withindex$Re.mosq.juv <- abs(as.numeric(Re_sim[,1]))
  traject.withindex$Re.mosq.adu <- abs(as.numeric(Re_sim[,2]))
  traject.withindex$Re.mosq.resjuv <- abs(as.numeric(Re_sim[,3]))
  traject.withindex$Re.mosq.resadu <- abs(as.numeric(Re_sim[,4]))
  traject.withindex$Re.juv.mosq <- abs(as.numeric(Re_sim[,5]))
  traject.withindex$Re.adu.mosq <- abs(as.numeric(Re_sim[,6]))
  traject.withindex$Re.resjuv.mosq <- abs(as.numeric(Re_sim[,7]))
  traject.withindex$Re.resadu.mosq <- abs(as.numeric(Re_sim[,8]))
  traject.withindex$Re.juv <- abs(as.numeric(Re_sim[,9]))
  traject.withindex$Re.adu <- abs(as.numeric(Re_sim[,10]))
  traject.withindex$Re.bb <- abs(as.numeric(Re_sim[,11]))
  traject.withindex$Re.resjuv <- abs(as.numeric(Re_sim[,12]))
  traject.withindex$Re.resadu <- abs(as.numeric(Re_sim[,13]))
  traject.withindex$Re.res <- abs(as.numeric(Re_sim[,14]))
  traject.withindex$Re <- abs(as.numeric(Re_sim[,15]))
  traject.withindex$Re_7d <- rollmeanr(traject.withindex$Re, 7, fill=NA)


  traject.withindex$prev.Juv <- as.numeric(prev[,1])
  traject.withindex$prev.Adu <- as.numeric(prev[,2])
  traject.withindex$prev.both <- as.numeric(prev[,3])
  traject.withindex$prev.mosq <- as.numeric(prev[,4])

  traject.withindex$sero.Juv <- as.numeric(sero[,1])
  traject.withindex$sero.Adu <- as.numeric(sero[,2])
  traject.withindex$sero.bothbb <- as.numeric(sero[,3])
  traject.withindex$sero.resjuv <- as.numeric(sero[,4])
  traject.withindex$sero.resadu <- as.numeric(sero[,5])
  traject.withindex$sero.bothres <- as.numeric(sero[,6])

  traject.withindex$dead.Juv <- as.numeric(dead[,1])
  traject.withindex$dead.Adu <- as.numeric(dead[,2])
  traject.withindex$dead.both <- as.numeric(dead[,3])
  traject.withindex$Nbb <- as.numeric(biting[,1])
  traject.withindex$Nres <- as.numeric(biting[,2])
  traject.withindex$prop.bites.bb.vs.res <- as.numeric(biting[,3])
  traject.withindex$prop.bites.bb.vs.total <- as.numeric(biting[,4])
  traject.withindex$res.vsbb <- as.numeric(biting[,5])

  simOutput[[run]] <- traject.withindex

}

saveRDS(simOutput, file=sprintf(paste(opt$output, '/simOutput.rds', sep='')))


# * Space: Calculate R0, Re, prevalence --------------------------------------------

simOutput <- list()

for (run in 1:length(listofresults)) { 

  print(run)
  R0_sim <- list()
  prop.S <- list()
  Re_sim <- list()

  modeloutput <- listofresults[[run]][['modeloutput']]

  modeloutput <- modeloutput %>%
    mutate(date = as.Date(time, origin="2015-12-31")) %>%
    mutate(month = as.numeric(format(date,'%m'))) %>%
    mutate(year = as.numeric(format(date,'%Y'))) %>%
    filter(month>=4 & month<=9)

  # aggregate over time
  setDT(modeloutput)
  modeloutput <- rollup(modeloutput, j=lapply(.SD, mean), by=c("node","year")) # takes a long time
  # this function adds a final row with the colMeans -> needs to be removed
  modeloutput <- modeloutput[1:(nrow(modeloutput)-1),]

  traject.withindex <- prep.modeloutput(modeloutput=modeloutput)

  for(row in 1:nrow(traject.withindex)) {

    R0_sim_row <- calculate.R0(output=traject.withindex[row,], posterior=listofresults[[run]][[2]],
                               transmissionProbMH=0.88, IIR=0.67, mortality.adu=0.0011, mortality.juv=0.0059, recovery=0.25,
                               version="multi-host")
    R0_sim <- rbind(R0_sim, R0_sim_row)

    prop.S_row <- calculate.propS(output=traject.withindex[row,])
    prop.S <- rbind(prop.S, prop.S_row)

    Re_sim_row <- calculate.Re(R0_sim = R0_sim[row,], prop.S = prop.S[row,])
    Re_sim <- rbind(Re_sim, Re_sim_row)
  }

  traject.withindex$R0 <- abs(as.numeric(R0_sim[,15]))
  traject.withindex$R0_7d <- rollmeanr(traject.withindex$R0, 7, fill=NA)

  traject.withindex$prop.sus.juv <- as.numeric(prop.S[,1])
  traject.withindex$prop.sus.adu <- as.numeric(prop.S[,2])
  traject.withindex$prop.sus.resjuv <- as.numeric(prop.S[,3])
  traject.withindex$prop.sus.resadu <- as.numeric(prop.S[,4])
  traject.withindex$prop.sus.mosq <- as.numeric(prop.S[,5])

  traject.withindex$Re <- abs(as.numeric(Re_sim[,15]))
  traject.withindex$Re_7d <- rollmeanr(traject.withindex$Re, 7, fill=NA)

  simOutput[[run]] <- traject.withindex

}

saveRDS(simOutput, file=sprintf(paste(opt$output, '/simOutput_space.rds', sep='')))

