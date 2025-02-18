# ABC-SMC DEFAULT PARAMETERS


# specified in userParameters
MODEL_DEFINITION_FILEPATH <- NA # Must be defined by the user
MODEL_LIST <- list() # Must be defined by the user
PRIOR_DIST <- list() # Must be defined by the user

MAX_NUMBER_OF_GEN <- 15
NB_ACC_PRTCL_PER_GEN <- 1000
VAR_PRTRBTN_KRNL <- NA
VAR_PRTRBTN_KRNL_MIN <- 0.01
SCALING_SD_PRTRBTN_KRNL <- 1.0
MODEL_JUMP_PROB <- 0.1
NB_THRESHOLD <- 1
NEW_THRESHOLD_QUANTILE <- 0.8
DISTANCE_THRESHOLD_MIN <- 0.01
MAXATTEMPTS <- 100000
ACCEPTANCE_RATE_MIN <- 0.01
USE_LHS_FOR_FIRST_ITER <- TRUE

ACCEPTED_PARTICLES_FILEPATH <- "res/csv/all_accepted_particles.csv"
LAST_ACCEPTED_PARTICLES_FILEPATH <- "res/csv/last_accepted_particles.csv"
ALL_PARTICLES_FILEPATH <- "res/csv/all_particles.csv"
THRESHOLDS_FILEPATH <- "res/csv/thresholds.csv"
LHS_USED_FOR_FIRST_GEN_FILEPATH <- "res/csv/lhs_for_first_gen.csv"

TMP_ACCEPTED_PARTICLES_FILEPATH <- "../Output/abc/tmp/tmp_current_gen_accepted_particles"
TMP_ALL_TESTED_PARTICLES_FILEPATH <- "../Output/abc/tmp/tmp_current_gen_all_tested_particles"
TMP_ARRAY_JOB_STD_OUT <- "../Output/abc/tmp/std_out/"
TMP_ARRAY_JOB_STD_ERR <- "../Output/abc/tmp/std_err/"

SSOBS_FILEPATH <- "ssobs.csv"

ON_CLUSTER <- FALSE
CLUSTER_TYPE <- "slurm" # "sge"
SLURM_SCRIPT_TEMPLATE <- ''
SGE_SCRIPT_TEMPLATE <- ''
BASH_SCRIPT_TEMPLATE <- 'Rscript abc-smc_modelselection_parallel_subjob.R --gen=%d --jobid=%d --abcparam=%s'
MAX_CONCURRENT_JOBS <- 100


# default parameters. not specified elsewhere
DISTNAMES <- paste0("dist", as.character(seq(1,NB_THRESHOLD,1)))
MODEL_NAMES <- c()
PARAM_NAMES <- c()
COLUMN_NAMES <- c()

EXPERIMENT_FOLDERPATH <- ""








