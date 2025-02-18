# ABC-SMC USER PARAMETERS

runID <- "id_20240923"
indexnr <- ".1"
scenarioName <- 1

print(scenarioName)
print(runID)
print(indexnr)


# select which model to use for identifiability
MODEL_DEFINITION_FILEPATH <- "abc_defineModel_id.R"
  
PRIOR_DIST <- list("baseline" = list(c("vertical_transm", "unif", 0, 0.3), c("transmissionProbHM", "unif", 0, 1),
                                       c("scalingParameter", "unif", 0.03, 2), c("disIndMortality", "unif", 0.1, 1),
                                       c("FOI.n", "unif", 0, 0.3), c("FOI.m", "unif", 0, 0.3), c("FOI.s", "unif", 0, 0.3)))

# ABC settings
MAX_NUMBER_OF_GEN <- 10
NB_ACC_PRTCL_PER_GEN <- 200
VAR_PRTRBTN_KRNL <- NA
VAR_PRTRBTN_KRNL_MIN <- 0.01
SCALING_SD_PRTRBTN_KRNL <- 0.5
MODEL_JUMP_PROB <- 0.1
NB_THRESHOLD <- 5
NEW_THRESHOLD_QUANTILE <- 0.7
DISTANCE_THRESHOLD_MIN <- 0.01
MAXATTEMPTS <- 1000000
ACCEPTANCE_RATE_MIN <- 0.0001
USE_LHS_FOR_FIRST_ITER <- TRUE
NB_ITERATIONS <- 3

# load summary statistics from simulated data & info to calculate patterns
load("../Data/abc_identifiability/summStat.Rda")

load("../Data/abc_summaryStatistics/PCRprev_likelihood.Rda")        # prevalence in region N/M/S across year
load("../Data/abc_summaryStatistics/sero_likelihood.Rda")           # observation process sero data
load("../Data/abc_summaryStatistics/dead.prev_likelihood.Rda")      # prevalence in region N/M/S across year

# select which ss_obs dataset to use
ss_obs <- ss_obs[[scenarioName]]                                    # data belonging to parameter set 1

# create output filepaths
ACCEPTED_PARTICLES_FILEPATH <- paste0("../Output/abc/", runID, "/all_accepted_particles_", scenarioName, indexnr, ".csv", sep="")
LAST_ACCEPTED_PARTICLES_FILEPATH <- paste0("../Output/abc/", runID, "/last_accepted_particles_", scenarioName, indexnr, ".csv", sep="")
ALL_PARTICLES_FILEPATH <- paste0("../Output/abc/", runID, "/all_particles_", scenarioName, indexnr, ".csv",sep="")
THRESHOLDS_FILEPATH <- paste0("../Output/abc/", runID, "/thresholds_", scenarioName, indexnr, ".csv",sep="")
LHS_USED_FOR_FIRST_GEN_FILEPATH <- paste0("../Output/abc/", runID, "/lhs_for_first_gen_", scenarioName, indexnr, ".csv", sep="")

TMP_ACCEPTED_PARTICLES_FILEPATH <- paste0("../Output/abc/", runID, "/tmp_current_gen_accepted_particles_", scenarioName, indexnr, sep="")
TMP_ALL_TESTED_PARTICLES_FILEPATH <- paste0("../Output/abc/", runID, "/tmp_current_gen_all_tested_particles_", scenarioName, indexnr, sep="")
TMP_ARRAY_JOB_STD_OUT <- paste0("../Output/abc/", runID, "/std_out_", scenarioName, indexnr, sep="")
TMP_ARRAY_JOB_STD_ERR <- paste0("../Output/abc/", runID, "/std_err_", scenarioName, indexnr, sep="")


# settings for HPC
ON_CLUSTER <- TRUE # if FALSE, will run on local machine launching concurrent processes (MAX_CONCURRENT_JOBS)
CLUSTER_TYPE <- "slurm"
MAX_CONCURRENT_JOBS <- 256


SLURM_SCRIPT_TEMPLATE <- '#!/bin/bash
#SBATCH --job-name=job-array   # nom du job
#SBATCH --partition=fat_rome
#SBATCH --ntasks=256
#SBATCH --time=08:00:00
#SBATCH --output=output_%%A_%%a.out
#SBATCH --error=error_output_%%A_%%a.out

module load 2023
module load R/4.3.2-gfbf-2023a

echo "%s %s %d"

for i in `seq 1 256` ; do
    Rscript abc_subjob.R --gen=%d --jobid=$i --abcparam=%s &
done
wait
'

