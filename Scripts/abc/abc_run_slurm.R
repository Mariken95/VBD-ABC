# Test of ABC-SMC with parallel implementation

source("abc_function_parallel.R")

# Test with a very simple toy model

# source("example_abc-smc_SIR_models_base.R")
# source("abc_defineModel_sc40.baseline.R")


# create a reference trajectory
#set.seed(123456789)
# params = c("mu" = 0.025, "beta" = 1.5, "gamma" = 0.5)
# sum_stat_obs = SIR_tau_leap(params)
# # print(sum_stat_obs)
# write.csv(sum_stat_obs, "ssobs.csv", row.names=FALSE, quote=FALSE)

args=(commandArgs(TRUE))
print(args)
if (length (args) == 0) {
  print("No arguments supplied.")
} else {
  for (i in 1:length(args)) {
    eval (parse (text = args[[i]] ))
  }
}

# versionArguments <- cbind(scenarioName, runID, indexnr)
# save(versionArguments, file="versionArguments.RData")

# .abc_user_param_file_path <- userParameterFile

if (!is.na(previous_gens)) {
  previous_gens <- read.csv(previous_gens)
  previous_epsilons <- read.csv(previous_epsilons)
}

# example with multiple distances / thresholds
# res = abcsmc(.abc_user_param_file_path="abc_userParameters.R", .verbose=TRUE)
res = abcsmc(.abc_user_param_file_path=userParameterFile, .previous_gens = previous_gens, .previous_epsilons = previous_epsilons, .verbose=TRUE)
all_accepted_particles = res$particles
all_thresholds = res$thresholds




