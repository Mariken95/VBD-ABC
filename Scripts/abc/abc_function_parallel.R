

#  _ _ _
# | (_) |__  _ __ __ _ _ __ _   _
# | | | '_ \| '__/ _` | '__| | | |
# | | | |_) | | | (_| | |  | |_| |
# |_|_|_.__/|_|  \__,_|_|   \__, |
#                           |___/

# lpackages <- list("lhs", "dplyr", "scales", "ggplot2", "GGally", "ggridges", "RColorBrewer", "progress", "processx", "optparse", "filelock")
# for (pack in lpackages) {
#     if(!require(pack,character.only = T)) {
#         install.packages(pack)
#         require(pack,character.only = T)
#     }
# }

# # Algorithm-required packages
# library(lhs,lib="/home/WUR/wit109/R/library")
# library(tidyverse,lib="/home/WUR/wit109/R/library")
# library(scales,lib="/home/WUR/wit109/R/library")
# library(progress,lib="/home/WUR/wit109/R/library") 
# library(processx,lib="/home/WUR/wit109/R/library") 
# library(optparse,lib="/home/WUR/wit109/R/library")
# library(filelock,lib="/home/WUR/wit109/R/library")
# 
# # Model-specific pacakges
# library(SimInf,lib="/home/WUR/wit109/R/library")
# library(philentropy,lib="/home/WUR/wit109/R/library")
# library(lubridate,lib="/home/WUR/wit109/R/library")
# library(data.table,lib="/home/WUR/wit109/R/library")

# Algorithm-required packages
.libPaths( c("/home/mdwit/R/library" , .libPaths() ) )

library(lhs,lib="/home/mdwit/R/library")
library(tidyverse,lib="/home/mdwit/R/library")
library(scales,lib="/home/mdwit/R/library")
library(progress,lib="/home/mdwit/R/library") 
library(processx,lib="/home/mdwit/R/library") 
library(optparse,lib="/home/mdwit/R/library")
library(filelock,lib="/home/mdwit/R/library")

# Model-specific pacakges
library(SimInf,lib="/home/mdwit/R/library")
library(philentropy,lib="/home/mdwit/R/library")
library(lubridate,lib="/home/mdwit/R/library")
library(data.table,lib="/home/mdwit/R/library")


source("abc_defaultParameters.R")


#        _
#   __ _| |__   ___   ___ _ __ ___   ___
#  / _` | '_ \ / __| / __| '_ ` _ \ / __|
# | (_| | |_) | (__  \__ \ | | | | | (__
#  \__,_|_.__/ \___| |___/_| |_| |_|\___|


createParticle = function(.gen, .rep, .lhs_for_first_iter, .previous_acc_particles, .empirical_sd, .use_lhs_for_first_iter, .prior_dist, .model_names, .model_jump_prob){
    the_particle = list()
    if (.gen == 1) {
        # sample from prior of get from lhs
        if (.use_lhs_for_first_iter) {
            the_particle = as.list(.lhs_for_first_iter[.rep,])
        } else {
            the_particle[["model"]] = sample(.model_names, 1)
            for (pp in .prior_dist[[the_particle[["model"]]]]) {
                if (pp[2] == "unif") {
                    val = runif(1, min = as.double(pp[3]), max = as.double(pp[4]))
                    the_particle[[pp[1]]] = val
                } else {
                    stop(p[2] + "- type of distribution not supported in the current version")
                }
            }
        }
    } else {
        # sample a model from previous population
        the_particle[["model"]] = sample(.previous_acc_particles$model, 1)
        # perturb the model
        if (rbinom(1,1,.model_jump_prob)) {
            the_particle[["model"]] = sample(unique(.previous_acc_particles$model), 1)
        }
        # select all particle corresponding to the chosen model, in the previous population
        previous_acc_particles_subpop = .previous_acc_particles[.previous_acc_particles$model == the_particle[["model"]],]
        # sample from previous population with associated weight
        samp_idx = sample(seq_len(nrow(previous_acc_particles_subpop)), 1, prob=previous_acc_particles_subpop$pWeight)
        sampled_particle = previous_acc_particles_subpop[samp_idx, ]
        the_particle = sampled_particle[c("model",sapply(.prior_dist[[the_particle[["model"]]]],"[[",1))]
        #print(the_particle)
        # perturb the particle
        for (pp in .prior_dist[[the_particle[["model"]]]]) {
            if (pp[2] == "unif") {
                while (1!=0) {
                	#print(.empirical_sd)
                    prtrbd_val = rnorm(1, mean = the_particle[[pp[1]]],
                                            sd = .empirical_sd[[the_particle[["model"]]]][[pp[1]]])
                    # check if new value contained in initial prior
                    if (dunif(prtrbd_val, min = as.double(pp[3]), max = as.double(pp[4])) > 0) {
                        the_particle[[pp[1]]] = prtrbd_val
                        break
                    }
                }
            } else {
                stop(p[2] + " - type of distribution not supported in the current version")
            }
        }
    }
    return(the_particle)
}


computeWeight = function(.proposed_particle, .previous_acc_particles, .empirical_sd, .prior_dist){
    model_used = .proposed_particle[["model"]]
    pWeight = 0
    numerator = 1
    denominator = .previous_acc_particles[.previous_acc_particles$model == model_used,"pWeight"]
    for (pp in .prior_dist[[model_used]]) {
        if (pp[2] == "unif") {
            numerator = numerator * dunif(.proposed_particle[[pp[1]]], min = as.double(pp[3]), max = as.double(pp[4]))
            denominator = denominator * dnorm(.proposed_particle[[pp[1]]], mean = .previous_acc_particles[.previous_acc_particles$model == model_used,pp[1]], sd = .empirical_sd[[model_used]][[pp[1]]])
        } else {
            stop(pp[2] + " - type of distribution not supported in the current version")
        }
    }
    pWeight = numerator / sum(denominator)
    return(pWeight)
}


defineNextThreshold = function(.acc_particles, .nb_threshold, .distnames, .new_threshold_quantile){
    the_epsilon = c()
    if (.nb_threshold > 1) {
        the_epsilon = unname(apply( .acc_particles[,.distnames], 2, quantile, probs = .new_threshold_quantile, na.rm = TRUE))
    } else {
        the_epsilon = c(unname(quantile(.acc_particles$dist, .new_threshold_quantile)))
    }
    return(the_epsilon)
}


setEmpiricalSd = function(.acc_particles, .model_names, .prior_dist, .var_prtrbtn_krnl, .scaling_sd_prtrbtn_krnl, .var_prtrbtn_krnl_min){
    the_empirical_sd = list()
    for (mm in .model_names) {
        the_empirical_sd[[mm]] = list()
        for (pp in .prior_dist[[mm]]) {
            the_empirical_sd[[mm]][[pp[1]]] = .var_prtrbtn_krnl
        }
    }
    #
    if (is.na(.var_prtrbtn_krnl)) {
        for (mm in .model_names) {
            for (pp in .prior_dist[[mm]]) {
                if (length(.acc_particles[.acc_particles$model == mm,pp[1]])>1) {
                    the_empirical_sd[[mm]][[pp[1]]] = sd(.acc_particles[.acc_particles$model == mm,pp[1]]) * .scaling_sd_prtrbtn_krnl[[pp[1]]]
                } else {
                    the_empirical_sd[[mm]][[pp[1]]] = .var_prtrbtn_krnl_min
                }
            }
        }
    }
    return(the_empirical_sd)
}


createLHSfromPrior = function(.model_names, .prior_dist, .nb_acc_prtcl_per_gen){
    all_lhs_for_first_iter = list()
    for (mm in .model_names) {
        tmp_lhs = as.data.frame(randomLHS(.nb_acc_prtcl_per_gen, length(.prior_dist[[mm]])))
        the_param_names = sapply(.prior_dist[[mm]], `[[`, 1)
        colnames(tmp_lhs) = the_param_names
        tmp_lhs["model"] = mm
        for (pp in .prior_dist[[mm]]) {
            tmp_lhs[,pp[1]] = rescale(tmp_lhs[,pp[1]], to = c(as.double(pp[3]),as.double(pp[4])), from = c(0.0,1.0))
        }
        all_lhs_for_first_iter[[mm]] = tmp_lhs
    }
    the_lhs_for_first_iter = bind_rows(all_lhs_for_first_iter)
    return(the_lhs_for_first_iter)
}


abcsmc = function(.abc_user_param_file_path, .previous_gens = NA, .previous_epsilons = NA, .verbose = FALSE){
    #
    source(.abc_user_param_file_path)
    source(MODEL_DEFINITION_FILEPATH)
    source("abc_updateParameters.R")

    # TODO : use a main folder with default name based on the datehoursminutesseconds if not specified
    for (file_path in c(TMP_ACCEPTED_PARTICLES_FILEPATH, TMP_ALL_TESTED_PARTICLES_FILEPATH, ACCEPTED_PARTICLES_FILEPATH, ALL_PARTICLES_FILEPATH, THRESHOLDS_FILEPATH, SSOBS_FILEPATH, LHS_USED_FOR_FIRST_GEN_FILEPATH, TMP_ARRAY_JOB_STD_OUT, TMP_ARRAY_JOB_STD_ERR)) {
        if (.verbose) {cat(paste0("Check filepath for : ", file_path, "\n"))}
        # Extract the directory part of the file path
        folder_path <- dirname(file_path)
        # Check if the folder path exists
        if (!dir.exists(folder_path)) {
          # Folder does not exist, create the folder
          if (dir.create(folder_path, recursive = TRUE)) {
            if (.verbose) {cat("Folder created successfully.\n")}
          } else {
            if (.verbose) {cat("Error: Unable to create folder.\n")}
          }
        } else {
          if (.verbose) {cat("Folder already exists.\n")}
        }
    }

    #
    gen = 1
    epsilon_not_improved = 0

    all_acc_particles = setNames(data.frame(matrix(ncol = length(COLUMN_NAMES), nrow = 0), stringsAsFactors=FALSE), COLUMN_NAMES)
    previous_acc_particles = setNames(data.frame(matrix(ncol = length(COLUMN_NAMES), nrow = 0), stringsAsFactors=FALSE), COLUMN_NAMES)
    current_acc_particles = setNames(data.frame(matrix(ncol = length(COLUMN_NAMES), nrow = 0), stringsAsFactors=FALSE), COLUMN_NAMES)

    epsilons = setNames(data.frame(matrix(ncol = length(DISTNAMES)+2, nrow = 0), stringsAsFactors=FALSE), c("gen", DISTNAMES, "acceptance_rate"))

    lhs_for_first_iter = c()
    if (USE_LHS_FOR_FIRST_ITER) {
        if (!all_prior_are_unif) {
            USE_LHS_FOR_FIRST_ITER <- FALSE
        } else {
            lhs_for_first_iter = createLHSfromPrior(MODEL_NAMES, PRIOR_DIST, NB_ACC_PRTCL_PER_GEN)
            write.csv(lhs_for_first_iter, LHS_USED_FOR_FIRST_GEN_FILEPATH, row.names=FALSE, quote=FALSE)
        }
    }

    # init with previous results if provided
    if (!is.null(dim(.previous_gens))) {
        # all_acc_particles
        all_acc_particles = .previous_gens
        epsilons = .previous_epsilons
        gen = max(all_acc_particles$gen)
        # previous_acc_particles : set based on last in all_acc_particles
        previous_acc_particles = as.data.frame(all_acc_particles[all_acc_particles$gen == gen,])
        # current_acc_particles  # nothing to do
        # epsilon
        epsilon = defineNextThreshold(previous_acc_particles, NB_THRESHOLD, DISTNAMES, NEW_THRESHOLD_QUANTILE)
        # calculate the empirical sd to be used in the perturbation kernel
        empirical_sd = setEmpiricalSd(previous_acc_particles, MODEL_NAMES, PRIOR_DIST, VAR_PRTRBTN_KRNL, SCALING_SD_PRTRBTN_KRNL, VAR_PRTRBTN_KRNL_MIN)
        # set the generation number
        gen = gen + 1
    }
    #
    while (gen <= MAX_NUMBER_OF_GEN) {
        # rep = 1
        nb_accepted = 0
        totattempts = 0
        if (.verbose) {
            cat("gen", gen, '\n')
            cat("threshold:", epsilon, '\n')
            cat("prtrbtn_krnl_sd:", unlist(empirical_sd), '\n')
        }
        #
        tmp_acc_prtcls = setNames(data.frame(matrix(ncol = length(COLUMN_NAMES), nrow = 0), stringsAsFactors=FALSE), COLUMN_NAMES)
        tmp_all_prtcls = setNames(data.frame(matrix(ncol = length(COLUMN_NAMES), nrow = 0), stringsAsFactors=FALSE), COLUMN_NAMES)
        write.csv(tmp_acc_prtcls, paste0(TMP_ACCEPTED_PARTICLES_FILEPATH,"_",gen,".csv"), row.names=FALSE, quote=FALSE)
        write.csv(tmp_all_prtcls, paste0(TMP_ALL_TESTED_PARTICLES_FILEPATH,"_",gen,".csv"), row.names=FALSE, quote=FALSE)
        #
        cluster_job_id = NA
        local_processes = c()
        local_processes_pids = c()
        if (ON_CLUSTER) {
            # submit a job array with as many job as N_acc_prtcl_before_next_gen (for gen == 1, need at one job per each particle)
            cluster_script = ""
            # Replace placeholders in the template with actual parameter values
            nbjobs = MAX_CONCURRENT_JOBS
            if (gen == 1) {nbjobs = min(N_acc_prtcl_before_next_gen, MAX_CONCURRENT_JOBS)}
            if (CLUSTER_TYPE == "slurm") {
                cluster_script <- sprintf(SLURM_SCRIPT_TEMPLATE, 1, nbjobs, MAX_CONCURRENT_JOBS, gen, .abc_user_param_file_path)
            } else if (CLUSTER_TYPE == "sge") {
                cluster_script <- sprintf(SGE_SCRIPT_TEMPLATE, 1, nbjobs, MAX_CONCURRENT_JOBS, gen, .abc_user_param_file_path)
            } else {
                stop("Cluster type not supported in the current version")
            }
            # Write the modified cluster job script to a file
            writeLines(cluster_script, "abc_smc_array_job.sh")
            # Submit the array job and capture the job ID
            if (CLUSTER_TYPE == "slurm") {
                cluster_job_id <- system('sbatch abc_smc_array_job.sh | awk \'{print $4}\'', intern = TRUE)
                print(paste("Submitted Slurm array job with ID:", cluster_job_id))
            } else if (CLUSTER_TYPE == "sge") {
                cluster_job_id <- system('qsub abc_smc_array_job.sh | awk -F "." \'{print $1}\' | awk \'{print $3}\'', intern = TRUE)
                print(paste("Submitted SGE array job with ID:", cluster_job_id))
            } else {
                stop("Cluster type not supported in the current version")
            }
        } else {
            # Parallel task on local machine
            for (task_index in 1:MAX_CONCURRENT_JOBS) {
                # Launch external script
                bash_script <- sprintf(BASH_SCRIPT_TEMPLATE, gen, task_index, .abc_user_param_file_path)
                writeLines(bash_script, paste0("abc_smc_task_", task_index, ".sh"))
                system(paste0("chmod 755 abc_smc_task_", task_index, ".sh"))
                command <- paste0("./abc_smc_task_", task_index, ".sh")
                local_process = process$new(command, stdout = "|", stderr = "|")
                local_processes <- c(local_processes, local_process)
                # Get the process id
                local_process_pid <- local_process$get_pid()
                local_processes_pids <- c(local_processes_pids, local_process_pid)
            }
            cat("local processes pids:", local_processes_pids, "\n")
        }
        #
        # Create a progress bar
        pb <- progress_bar$new(format = "gen :gen [:bar] :percent (ar: :accrate | :nbattempt) | eta: :eta (:elapsed)", clear = FALSE, total = N_acc_prtcl_before_next_gen)
        #
        current_iter_broke = FALSE
        while (nb_accepted < N_acc_prtcl_before_next_gen) { # repeat until N particles accepted
            tmp_accepted_particles = read.csv(paste0(TMP_ACCEPTED_PARTICLES_FILEPATH,"_",gen,".csv"))
            tmp_all_tested_particles = read.csv(paste0(TMP_ALL_TESTED_PARTICLES_FILEPATH,"_",gen,".csv"))
            nb_accepted = nrow(tmp_accepted_particles)
            totattempts = nrow(tmp_all_tested_particles)
            current_acc_rate = nb_accepted/totattempts
            # Increment the progress bar
            pb$update(min(nb_accepted/N_acc_prtcl_before_next_gen, 1.00), tokens = list(gen = gen, accrate = format(round(current_acc_rate,digits=3),nsmall=3), nbattempt = totattempts))
            #
            if (totattempts > MAXATTEMPTS) {
                message('\n', "The maximum number of attempts to accept the N particles has been reached!", '\n')
                current_iter_broke = TRUE
                break
            }
            if ((totattempts >= (NB_ACC_PRTCL_PER_GEN*length(MODEL_NAMES))) & (current_acc_rate < ACCEPTANCE_RATE_MIN)) {
                message('\n', "The acceptance rate has become too low (< specified acceptance_rate_min), the algorithm stops!")
                current_iter_broke = TRUE
                break
            }
        }
        # cancel subjob on cluster
        if (gen == 1) {
            # nothing to do
        } else {
            if (ON_CLUSTER) {
                # Cancel the array job
                if (CLUSTER_TYPE == "slurm") {
                    cancel_command <- paste("scancel", cluster_job_id)
                    system(cancel_command, intern = TRUE)
                } else if (CLUSTER_TYPE == "sge") {
                    cancel_command <- paste("qdel", cluster_job_id)
                    system(cancel_command, intern = TRUE)
                } else {
                    stop("Cluster type not supported in the current version")
                }
            } else {
                # kill local process if needed
                for (local_process in local_processes) {
                    if (local_process$is_alive()) {
                        local_process$kill()
                    }
                }
            }
        }
        # Close the progress bar
        # pb$terminate()
        # cat('\n')
        #
        if (current_iter_broke == TRUE) {break}
        #
        current_acc_particles = read.csv(paste0(TMP_ACCEPTED_PARTICLES_FILEPATH,"_",gen,".csv"))
        current_acc_particles = current_acc_particles[1:min(nrow(current_acc_particles),N_acc_prtcl_before_next_gen),] # keep only the number of particle needed # TODO : imporve comment
        # normalise the weights
        for (mm in MODEL_NAMES) {
            if (mm %in% unique(current_acc_particles$model)) {
        	   current_acc_particles[current_acc_particles$model == mm,"pWeight"]= current_acc_particles[current_acc_particles$model == mm,"pWeight"]/ sum(current_acc_particles[current_acc_particles$model == mm,"pWeight"])
            }
        }
        write.csv(current_acc_particles, LAST_ACCEPTED_PARTICLES_FILEPATH, row.names=FALSE, quote=FALSE)
        # switch current and previous particles
        previous_acc_particles = current_acc_particles
        # save current particles
        all_acc_particles = rbind(all_acc_particles, current_acc_particles)
        # clear current particles table for next generation
        current_acc_particles = current_acc_particles[c(),]
        #
        if (gen > 1) {
            # save current epsilon if generation completed
            epsilons[nrow(epsilons) + 1,] = c(gen, epsilon, current_acc_rate)
        }
        # save results
        write.csv(all_acc_particles, ACCEPTED_PARTICLES_FILEPATH, row.names=FALSE, quote=FALSE)
        write.csv(epsilons, THRESHOLDS_FILEPATH, row.names=FALSE, quote=FALSE)
        # STOPPING RULES
        # Stop the algo if there is no significant difference between two successive posterior distributions
        # if (gen > 1) {
        #     # To be implemented
        # }
        # Stop the algo if epsilon falls below the predetermined min threshold
        if ((gen > 1) & (all(epsilon <= DISTANCE_THRESHOLD_MIN))) {
            message("The distance threshold(s) (epsilon(s)) fall(s) below the predetermined min value!")
            print(epsilon)
            break
        }
        # calculate the empirical sd to be used in the perturbation kernel
        empirical_sd = setEmpiricalSd(previous_acc_particles, MODEL_NAMES, PRIOR_DIST, VAR_PRTRBTN_KRNL, SCALING_SD_PRTRBTN_KRNL, VAR_PRTRBTN_KRNL_MIN)
        # define the next threshold
        next_epsilon = defineNextThreshold(previous_acc_particles, NB_THRESHOLD, DISTNAMES, NEW_THRESHOLD_QUANTILE)
        if ((gen > 1) & (all(next_epsilon == epsilon))) {
            epsilon_not_improved = epsilon_not_improved + 1
            if (epsilon_not_improved >= 2) {
                message("The distance threshold(s) (epsilon(s)) has not been improved over the last two iterations!")
                print(epsilons)
                break
            }
        }
        epsilon = next_epsilon
        #
        gen = gen + 1
        #
        # if (.verbose) {cat('\n-\n')}
        if (.verbose) {cat('-\n')}
    }
    if (.verbose) {
        cat("all thresholds:", '\n')
        print(epsilons)
        print(summary(previous_acc_particles[,c("model", PARAM_NAMES)]))
    }
    # cleaning
    if (ON_CLUSTER) {
        unlink("abc_smc_array_job.sh")
    } else {
        for (task_index in 1:MAX_CONCURRENT_JOBS) {
            unlink(paste0("abc_smc_task_", task_index, ".sh"))
        }
    }
    # rm tmp ?
    unlink("tmp", recursive = TRUE)
    #
    return(list("particles" = all_acc_particles, "thresholds" = epsilons))
}


#        _       _
#  _ __ | | ___ | |_
# | '_ \| |/ _ \| __|
# | |_) | | (_) | |_
# | .__/|_|\___/ \__|
# |_|


plot_abcsmc_res = function(data, prior, filename = "pairplot.pdf", figtitle = "ABC smc results", colorpal = "GnBu", iter = NA){
    # message("The plot may fail if, for one or more iterations, only one particle has been selected for one of the models.")
    tmp = data
    colorshift = 2  # because the first colours is often too light
    nb.cols = max(tmp$gen)+1+colorshift
    mycolors = colorRampPalette(brewer.pal(8, colorpal))(nb.cols)
    mycolors = mycolors[seq(-1,-colorshift,-1)]

    if (all(!is.na(iter))) {
        tmp = data[data$gen %in% iter,]
        mycolors = mycolors[iter]
    }
    tmp$gen = as.factor(tmp$gen)
    pdf(filename, width=9, height=9)
    # param_names = unique(Reduce(c, sapply(prior, function(x) sapply(x, `[[`, 1))))
    # pairplot = ggpairs(tmp[c(param_names, "gen")], aes(fill = gen, colour = gen, alpha=0.8),
    #     columns = param_names, lower = list(continuous = wrap("points", alpha = 0.5, size=1.5)),
    #     title = figtitle
    #  ) + scale_fill_manual(values = mycolors) + scale_colour_manual(values = mycolors) + theme_dark()
    # print(pairplot)
    if (length(unique(tmp$model)) >= 1) {
        if (length(unique(tmp$model)) > 1) {
            tmpm = tmp[,c("gen","model")]
            modelcolors = colorRampPalette(brewer.pal(8, colorpal))(length(unique(tmpm$model))+1)
            modelcolors = modelcolors[-1]
            stackedbarplot = ggplot(data, aes(x = factor(gen), fill = model)) +
              geom_bar(position = "fill") +
              labs(x = "gen", y = "Count", fill = "model") +
              ggtitle(paste(figtitle, 'model selection', sep=" - ")) +
              scale_fill_manual(values = modelcolors) +
              theme_minimal() +
              theme(legend.position="bottom", panel.grid.minor = element_blank(), panel.grid.major = element_blank(), panel.grid.major.y = element_line(linewidth=.1, color="black"))
            print(stackedbarplot)
        }
        for (mm in sort(unique(tmp$model))) {
            param_names = sapply(prior[[mm]],"[[",1)
            tmpm = tmp[tmp$model == mm,]
            pairplot = ggpairs(tmpm[c(param_names, "gen")], upper = "blank", aes(fill = gen, colour = gen, alpha=0.8),
                columns = param_names, lower = list(continuous = wrap("points", alpha = 0.5, size=1.5)),
                title = paste(figtitle, mm, sep=" - ")
             ) + scale_fill_manual(values = mycolors) + scale_colour_manual(values = mycolors) + theme_minimal() #theme_dark()
            print(pairplot)
        }
    }
    dev.off()
}


plot_densityridges = function(data, prior, filename = "densityridges.pdf", figtitle = "", colorpal = "Blues"){
    tmp = data
    colorshift = 2  # because the first colours is often too light
    nb.cols = max(tmp$gen)+1+colorshift
    mycolors = colorRampPalette(brewer.pal(8, colorpal))(nb.cols)
    mycolors = mycolors[seq(-1,-colorshift,-1)]
    tmp$gen = as.factor(tmp$gen)
    pdf(filename, width=6, height=8)
    for (mm in sort(unique(tmp$model))) {
        for (p in prior[[mm]]) {
            subtmp = tmp[tmp$model == mm,c('gen', p[1])]
            names(subtmp)[names(subtmp) == p[1]] = 'param'
            sub_plot = ggplot(subtmp, aes(x=param, y=gen, fill=gen, color=gen))+
                          geom_density_ridges(aes(), scale = 2.0, alpha = .6, bandwidth = (as.double(p[4]) - as.double(p[3]))/100) + theme_ridges() +
                          scale_y_discrete(expand = expansion(add = c(0.3, 2.5))) +
                          scale_x_continuous(limits = c(as.double(p[3]), as.double(p[4])), expand = c(0.01, 0)) +
                          scale_fill_manual(values = mycolors) +
                          scale_colour_manual(values = mycolors) +
                          labs(x=p[1], y="iteration")+
                          # ggtitle(paste("Successive posterior distributions for parameter",p[1], sep=' '))
                          theme_minimal() + theme(legend.position = "none")
            print(sub_plot)
        }
    }
    dev.off()
}


plot_thresholds = function(data, nb_threshold = 1, filename = "thresholds.pdf", figtitle = "", colorpal = "Blues"){
    tmp = data
    colorshift = 2  # because the first colours is often too light
    nb.cols = max(tmp$gen)+1+colorshift
    mycolors = colorRampPalette(brewer.pal(8, colorpal))(nb.cols)
    mycolors = mycolors[seq(-1,-colorshift,-1)]
    tmp$gen = as.factor(tmp$gen)
    pdf(filename, width=8, height=5)
    for (dd in paste0("dist", as.character(seq(1,nb_threshold,1)))) {
        subtmp = tmp[,c('gen', dd)]
        names(subtmp)[names(subtmp) == dd] = 'dist'
        sub_plot = ggplot(subtmp, aes(x=gen, y=dist, fill=gen, color=gen))+
                      geom_bar(stat = "identity") +
                      scale_fill_manual(values = mycolors) +
                      scale_colour_manual(values = mycolors) +
                      labs(y=paste0("threshold for ", dd), x="iteration ID") +
                      # ggtitle(paste("Successive posterior distributions for parameter",p[1], sep=' '))
                      theme_minimal() + theme(legend.position = "none")
        print(sub_plot)
    }
    dev.off()
}

###############################################################################
