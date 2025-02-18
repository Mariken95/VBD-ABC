
DISTNAMES <- paste0("dist", as.character(seq(1,NB_THRESHOLD,1)))
MODEL_NAMES <- names(MODEL_LIST)
PARAM_NAMES <- unique(Reduce(c, sapply(PRIOR_DIST, function(x) sapply(x, `[[`, 1))))
COLUMN_NAMES <- c("gen", "model", PARAM_NAMES, "pWeight", DISTNAMES)

# define the total number of particles to accept before next gen, that will be
# used as a upper limit for the number of simulation to run (avoid a while loop
# without stopping criterion)
N_acc_prtcl_before_next_gen <- NB_ACC_PRTCL_PER_GEN*length(MODEL_NAMES)

# define the vector of thresholds
epsilon = c()

# define the vector of sd to use in the perturbation kernel
empirical_sd = list()
for (mm in MODEL_NAMES) {
    empirical_sd[[mm]] = list()
    for (pp in PRIOR_DIST[[mm]]) {
        empirical_sd[[mm]][[pp[1]]] = VAR_PRTRBTN_KRNL
    }
}

if (length(SCALING_SD_PRTRBTN_KRNL) == 1) {
    scaling_sd_prtrbtn_krnl_uniqval = SCALING_SD_PRTRBTN_KRNL
    SCALING_SD_PRTRBTN_KRNL = list()
    for (pp in PARAM_NAMES) {
        SCALING_SD_PRTRBTN_KRNL[[pp]] = scaling_sd_prtrbtn_krnl_uniqval
    }
}

all_prior_are_unif = all(unique(Reduce(c, sapply(PRIOR_DIST, function(x) sapply(x, `[[`, 2)))) == "unif")
