# subjob.R

# Rscript subjob.R --gen=1 --jobid=42 --abcparam=abc_user_parameters.R


#  _                 _          _
# | | ___   __ _  __| |   _ __ | | ____ _ ___
# | |/ _ \ / _` |/ _` |  | '_ \| |/ / _` / __|
# | | (_) | (_| | (_| |  | |_) |   < (_| \__ \
# |_|\___/ \__,_|\__,_|  | .__/|_|\_\__, |___/
#                        |_|        |___/

source("abc_function_parallel.R")
source("abc_defaultParameters.R")


#                     _    _ _
#   ___ _ __ ___   __| |  | (_)_ __   ___     __ _ _ __ __ _ ___
#  / __| '_ ` _ \ / _` |  | | | '_ \ / _ \   / _` | '__/ _` / __|
# | (__| | | | | | (_| |  | | | | | |  __/  | (_| | | | (_| \__ \
#  \___|_| |_| |_|\__,_|  |_|_|_| |_|\___|   \__,_|_|  \__, |___/
#                                                      |___/

# Define command line options
option_list <- list(
  make_option(c("-a", "--abcparam"),
              type="character",
              default=NULL,
              help="abc user param filepath : an R file containing the ABC-SMC parameters defined by the user"),
  make_option(c("-g", "--gen"),
              type="integer",
              default=NULL,
              help="id of the current gen (integer)"),
  make_option(c("-j", "--jobid"),
              type="integer",
              default=NULL,
              help="id of the array job (integer)")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# id_current_gen
if (!is.null(opt$gen)) {
  id_current_gen <- opt$gen
  print(paste("the current gen:", id_current_gen))
}

# job_id
if (!is.null(opt$jobid)) {
  job_id <- opt$jobid
  print(paste("the job id:", job_id))
}

if (!is.null(opt$abcparam)) {
  abc_user_param_file_path <- opt$abcparam
  print(paste("abc user param file path :", abc_user_param_file_path))
}



#           _                      _       _     _
#  ___  ___| |_   __   ____ _ _ __(_) __ _| |__ | | ___  ___
# / __|/ _ \ __|  \ \ / / _` | '__| |/ _` | '_ \| |/ _ \/ __|
# \__ \  __/ |_    \ V / (_| | |  | | (_| | |_) | |  __/\__ \
# |___/\___|\__|    \_/ \__,_|_|  |_|\__,_|_.__/|_|\___||___/

# load arguments from python script
# load(file="versionArguments.RData")
# scenarioName <- versionArguments[1]
# runID <- versionArguments[2]
# indexnr <- versionArguments[3]

# source("abc_userParameters.R")
source(abc_user_param_file_path)
source(MODEL_DEFINITION_FILEPATH)
source("abc_updateParameters.R")


# define the dataframe to store previous accepted particles
previous_acc_particles = setNames(data.frame(matrix(ncol = length(COLUMN_NAMES), nrow = 0), stringsAsFactors=FALSE), COLUMN_NAMES)

# set some variables depending on the generation id
current_gen_acc_prtcls = 1
if (id_current_gen == 1) {
  current_gen_acc_prtcls = 1 + (job_id-1)*(N_acc_prtcl_before_next_gen%/%min(MAX_CONCURRENT_JOBS,N_acc_prtcl_before_next_gen))
  if (job_id < min(MAX_CONCURRENT_JOBS,N_acc_prtcl_before_next_gen)) {
    N_acc_prtcl_before_next_gen = current_gen_acc_prtcls + (N_acc_prtcl_before_next_gen%/%min(MAX_CONCURRENT_JOBS,N_acc_prtcl_before_next_gen)) - 1
  }
} else { # id_current_gen > 1
  # get previous accepted particles
  previous_gen <- read.csv(ACCEPTED_PARTICLES_FILEPATH)
  previous_acc_particles = as.data.frame(previous_gen[previous_gen$gen == id_current_gen-1,])
  # calculate the empirical sd to be used in the perturbation kernel
  empirical_sd = setEmpiricalSd(previous_acc_particles, MODEL_NAMES, PRIOR_DIST, VAR_PRTRBTN_KRNL, SCALING_SD_PRTRBTN_KRNL, VAR_PRTRBTN_KRNL_MIN)
  # define the current threshold
  epsilon = defineNextThreshold(previous_acc_particles, NB_THRESHOLD, DISTNAMES, NEW_THRESHOLD_QUANTILE)
}


#  _                 _        _       _
# | | ___   __ _  __| |    __| | __ _| |_ __ _
# | |/ _ \ / _` |/ _` |   / _` |/ _` | __/ _` |
# | | (_) | (_| | (_| |  | (_| | (_| | || (_| |
# |_|\___/ \__,_|\__,_|   \__,_|\__,_|\__\__,_|

# load the required data
# ss_obs <- read.csv(SSOBS_FILEPATH)
lhs_first_gen <- c()
if (USE_LHS_FOR_FIRST_ITER) {
  if (!all_prior_are_unif) {
    USE_LHS_FOR_FIRST_ITER <- FALSE
    } else {
      lhs_first_gen <- read.csv(LHS_USED_FOR_FIRST_GEN_FILEPATH)
    }
}


#  _            _                       _   _      _
# | |_ ___  ___| |_    _ __   __ _ _ __| |_(_) ___| | ___  ___
# | __/ _ \/ __| __|  | '_ \ / _` | '__| __| |/ __| |/ _ \/ __|
# | ||  __/\__ \ |_   | |_) | (_| | |  | |_| | (__| |  __/\__ \
#  \__\___||___/\__|  | .__/ \__,_|_|   \__|_|\___|_|\___||___/
#                     |_|

while (current_gen_acc_prtcls <= N_acc_prtcl_before_next_gen) {
  # create the particle
  proposed_particle = createParticle(id_current_gen, current_gen_acc_prtcls, lhs_first_gen, previous_acc_particles, empirical_sd, USE_LHS_FOR_FIRST_ITER, PRIOR_DIST, MODEL_NAMES, MODEL_JUMP_PROB)
  # simulate and compute the distance
  # dist = MODEL_LIST[[proposed_particle[["model"]]]](proposed_particle, ss_obs) # 1 iteration per particle
  dist = lapply(seq_len(NB_ITERATIONS), function(x) MODEL_LIST[[proposed_particle[["model"]]]](proposed_particle, ss_obs)) # multiple iterations per particle
  dist <- do.call("rbind", dist)
  print(dist)
  dist <- colMeans(dist)
  
  # compute the weight
  pWeight = NA
  if (id_current_gen == 1) {
      # compute the weight
      pWeight = 1
  } else {
      # accept particle if d <= epsilon, and so compute the weight
      if (all(dist <= epsilon)) {
          # compute the weight
          pWeight = computeWeight(proposed_particle, previous_acc_particles, empirical_sd, PRIOR_DIST)
      }
  }

  # save particle in the shared table
  if ((id_current_gen == 1) || (all(dist <= epsilon))) {
    # build the new row
    new.row = c(list(gen=id_current_gen, pWeight=pWeight), proposed_particle, setNames(as.list(dist), DISTNAMES))
    new.row = data.frame(new.row)
    missing_columns = setdiff(COLUMN_NAMES, colnames(new.row)) # add the missing columns to the new row with empty values (or NA)
    new.row[missing_columns] = NA  # You can also define other default values if required
    new.row = new.row[COLUMN_NAMES] # sort the columns to keep them in the right order
    # put lock on resfile.csv.lck
    lck = lock(paste0(TMP_ACCEPTED_PARTICLES_FILEPATH,"_",id_current_gen,".csv.lck"))
    # Write the new line to the CSV file without reading it first
    write.table(new.row, file = paste0(TMP_ACCEPTED_PARTICLES_FILEPATH,"_",id_current_gen,".csv"), sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
    # remove lock on resfile.csv.lck
    unlock(lck)
    # increment the number of accepted particle
    current_gen_acc_prtcls = current_gen_acc_prtcls + 1
  }

  # in any case, add the tested particle in a file
  # build the new row
  new.row = c(list(gen=id_current_gen, pWeight=pWeight), proposed_particle, setNames(as.list(dist), DISTNAMES))
  new.row = data.frame(new.row)
  missing_columns = setdiff(COLUMN_NAMES, colnames(new.row)) # add the missing columns to the new row with empty values (or NA)
  new.row[missing_columns] = NA  # You can also define other default values if required
  new.row = new.row[COLUMN_NAMES] # sort the columns to keep them in the right order
  # put lock on resfile.csv.lck
  lck = lock(paste0(TMP_ALL_TESTED_PARTICLES_FILEPATH,"_",id_current_gen,".csv.lck"))
  # Write the new line to the CSV file without reading it first
  write.table(new.row, file = paste0(TMP_ALL_TESTED_PARTICLES_FILEPATH,"_",id_current_gen,".csv"), sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
  # remove lock on resfile.csv.lck
  unlock(lck)

}
