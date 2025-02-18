

# General -----------------------------------------------------------------

combine_rds_files <- function(folder_paths, file_name) {
  combined_data <- list()
      # browser()

  for (folder in folder_paths) {
    file_path <- file.path(folder, file_name)
    
    data <- readRDS(file_path)
    
    # Assuming the loaded data is a data frame and is named `data`
    if (is.null(combined_data)) {
      combined_data <- data
    } else {
      combined_data <- c(combined_data, data)
    }
  } 
  
  return(combined_data)
}


# Summary statistics ------------------------------------------------------

get_ss.sim.output <- function(summ.stat.dimensions, ss_sim_posterior, patternnr) {
  
  simulations <- data.frame(matrix(nrow=summ.stat.dimensions, ncol=0))
  
  for (ii in 1:length(ss_sim_posterior)) {
    simulations <- cbind(simulations, ss_sim_posterior[[ii]][[patternnr]][3])
  }
  # add year and region column
  simulations <- cbind(simulations, ss_sim_posterior[[1]][[patternnr]][1:2])
  
  # browser()
  colnames(simulations) <- c(paste("sim", 1:length(ss_sim_posterior), sep="_"),"region", "year")
  
  simulations <- simulations %>%
    rowwise(.) %>%
    mutate(mean = mean(c_across(c(all_of(grep("sim", names(.), fixed=TRUE)))))) %>%
    mutate(lb = quantile(c_across(c(all_of(grep("sim", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub = quantile(c_across(c(all_of(grep("sim", names(.), fixed=TRUE)))), c(0.025, 0.975))[2])
  
  return(simulations)
}


# R0 ----------------------------------------------------------------------

calculate.R0 <- function(output, posterior, 
                         transmissionProbMH, IIR, mortality.adu, mortality.juv, recovery, 
                         version="multi-host") {
  
  if (version=="single-host") {
    # browser()
    Njuv <- output$Shj + output$Ehj + output$Ihj + output$Rhj
    Nadu <- output$Sha + output$Eha + output$Iha + output$Rha
    Nmos <- output$Sm + output$Em + output$Im 
    
    mosq.juv <- (output$biting * output$prop.bites * transmissionProbMH * output$incubationRateM * Njuv) / 
      ((output$incubationRateM + output$mumos + output$diapause) * (output$mumos + output$diapause) * (Njuv + Nadu))
    
    mosq.adu <- (output$biting * output$prop.bites * transmissionProbMH * output$incubationRateM * Nadu) / 
      ((output$incubationRateM + output$mumos + output$diapause) * (output$mumos + output$diapause) * (Njuv + Nadu))
    
    juv.mosq <- (output$biting * output$prop.bites * posterior$transmissionProbHM * IIR * Njuv * Nmos) /
      ((IIR + mortality.juv) * (mortality.juv + recovery + posterior$disIndMortality) * (Njuv + Nadu) * (Njuv + Nadu))
    
    adu.mosq <- (output$biting * output$prop.bites * posterior$transmissionProbHM * IIR * Nadu * Nmos) /
      ((IIR + mortality.adu) * (mortality.adu + recovery + posterior$disIndMortality) * (Njuv + Nadu) * (Njuv + Nadu))
    
    NGM <- matrix(data=c(0, juv.mosq, adu.mosq, 
                         mosq.juv, 0, 0, 
                         mosq.adu, 0, 0), 
                  nrow=3, ncol=3)
  }
  
  else if (version=="multi-host") {
    # browser()
    
    Njuv <- output$Shj + output$Ehj + output$Ihj + output$Rhj
    Nadu <- output$Sha + output$Eha + output$Iha + output$Rha
    Nresjuv <- output$Srj + output$Erj + output$Irj + output$Rrj
    Nresadu <- output$Sra + output$Era + output$Ira + output$Rra
    Nmos <- output$Sm + output$Em + output$Im 
    
    R0.mosq.juv <- ifelse(Njuv>0,
                          (output$incubationRateM * output$biting * output$prop.bites * transmissionProbMH * Njuv) /
                            ((output$incubationRateM + output$mumos + output$diapause) * (output$mumos + output$diapause) * ((Nresjuv + Nresadu) * posterior$bitingPref + Njuv + Nadu)),
                          0)
    
    R0.mosq.adu <- ifelse(Nadu>0,
                          (output$incubationRateM * output$biting * output$prop.bites * transmissionProbMH * Nadu) / 
                            ((output$incubationRateM + output$mumos + output$diapause) * (output$mumos + output$diapause) * ((Nresjuv + Nresadu) * posterior$bitingPref + Njuv + Nadu)),
                          0)
    
    R0.mosq.resjuv <- ifelse(Nresjuv>0,
                             (output$incubationRateM * output$biting * output$prop.bites * transmissionProbMH *  (Nresjuv * posterior$bitingPref)) / 
                               ((output$incubationRateM + output$mumos + output$diapause) * (output$mumos + output$diapause) * ((Nresjuv + Nresadu) * posterior$bitingPref + Njuv + Nadu)),
                             0)
    
    R0.mosq.resadu <- ifelse(Nresadu>0,
                             (output$incubationRateM * output$biting * output$prop.bites * transmissionProbMH *  (Nresadu * posterior$bitingPref)) / 
                               ((output$incubationRateM + output$mumos + output$diapause) * (output$mumos + output$diapause) * ((Nresjuv + Nresadu) * posterior$bitingPref + Njuv + Nadu)),
                             0)

    R0.juv.mosq <- ifelse(Njuv>0, 
                          (IIR * output$biting * output$prop.bites * posterior$transmissionProbHM * Njuv * Nmos) /
                            ((IIR + mortality.juv) * (mortality.juv + recovery + posterior$disIndMortality) * ((Nresjuv + Nresadu) * posterior$bitingPref + Njuv + Nadu) * Njuv),
                          0)
    
    R0.adu.mosq <- ifelse(Nadu>0,
                          (IIR * output$biting * output$prop.bites * posterior$transmissionProbHM * Nadu * Nmos) /
                            ((IIR + mortality.adu) * (mortality.adu + recovery + posterior$disIndMortality) * ((Nresjuv + Nresadu) * posterior$bitingPref + Njuv + Nadu) * Nadu),
                          0)
    
    R0.resjuv.mosq <- ifelse(Nresjuv>0,
                             (IIR * output$biting * output$prop.bites * posterior$transmissionProbHM * (Nresjuv * posterior$bitingPref) * Nmos) /
                               ((IIR + mortality.juv) * (mortality.juv + recovery) * ((Nresjuv + Nresadu) * posterior$bitingPref + Njuv + Nadu) * Nresjuv),
                             0)
    
    R0.resadu.mosq <- ifelse(Nresadu>0,
                             (IIR * output$biting * output$prop.bites * posterior$transmissionProbHM * (Nresadu * posterior$bitingPref) * Nmos) /
                               ((IIR + (1/posterior$lifespanR)) * ((1/posterior$lifespanR) + recovery) * ((Nresjuv + Nresadu) * posterior$bitingPref + Njuv + Nadu) * Nresadu),
                             0)

        NGM <- matrix(data=c(0, R0.juv.mosq, R0.adu.mosq, R0.resjuv.mosq, R0.resadu.mosq,
                            R0.mosq.juv, 0, 0, 0, 0,
                            R0.mosq.adu, 0, 0, 0, 0,
                            R0.mosq.resjuv, 0, 0, 0, 0,
                            R0.mosq.resadu, 0, 0, 0, 0), 
                  nrow=5, ncol=5)
  }
  
  eigenvalues <- eigen(NGM)
  R0 <- max(abs(eigenvalues$values))
  print(paste("R0 from NGM:", R0))
  
  R0.juv <- sqrt(R0.mosq.juv * R0.juv.mosq)
  R0.adu <- sqrt(R0.mosq.adu * R0.adu.mosq)
  R0.bb <- sqrt(R0.mosq.juv * R0.juv.mosq + R0.mosq.adu * R0.adu.mosq) 
  R0.resjuv <- sqrt(R0.mosq.resjuv * R0.resjuv.mosq)
  R0.resadu <- sqrt(R0.mosq.resadu * R0.resadu.mosq)
  R0.res <- sqrt(R0.mosq.resjuv * R0.resjuv.mosq + R0.mosq.resadu * R0.resadu.mosq)
  
  Routput <- cbind(R0.mosq.juv, R0.mosq.adu, R0.mosq.resjuv, R0.mosq.resadu, R0.juv.mosq, R0.adu.mosq, R0.resjuv.mosq, R0.resadu.mosq,
                   R0.juv, R0.adu, R0.bb, R0.resjuv, R0.resadu, R0.res, R0)
  return(Routput)
  
}

calculate.propS <- function(output) {
  
  # browser()
  Njuv <- output$Shj + output$Ehj + output$Ihj + output$Rhj
  Nadu <- output$Sha + output$Eha + output$Iha + output$Rha
  Nresjuv <- output$Srj + output$Erj + output$Irj + output$Rrj
  Nresadu <- output$Sra + output$Era + output$Ira + output$Rra
  Nmos <- output$Sm + output$Em + output$Im 
  
  prop.sus.juv <- ifelse(Njuv>0, (output$Shj / Njuv ), 0)
  prop.sus.adu <- ifelse(Nadu>0, (output$Sha / Nadu), 0)
  prop.sus.resjuv <- ifelse(Nresjuv>0, (output$Srj / Nresjuv), 0)
  prop.sus.resadu <- ifelse(Nresadu>0, (output$Sra / Nresadu), 0)
  prop.sus.mosq <- ifelse(Nmos>0, (output$Sm / Nmos), 0)
  
  prop.sus <- cbind(prop.sus.juv, prop.sus.adu, prop.sus.resjuv, prop.sus.resadu, prop.sus.mosq)
  
  print(prop.sus)
  return(prop.sus)
  
}

calculate.Re <- function(R0_sim, prop.S) {
  # browser()
  Re.mosq.juv <- R0_sim[[1]] * prop.S[[1]]
  Re.mosq.adu <- R0_sim[[2]] * prop.S[[2]]
  Re.mosq.resjuv <- R0_sim[[3]] * prop.S[[3]]
  Re.mosq.resadu <- R0_sim[[4]] * prop.S[[4]]

  Re.juv.mosq <- R0_sim[[5]] * prop.S[[5]]
  Re.adu.mosq <- R0_sim[[6]] * prop.S[[5]]
  Re.resjuv.mosq <- R0_sim[[7]] * prop.S[[5]]
  Re.resadu.mosq <- R0_sim[[8]] * prop.S[[5]]
  
  NGM <- matrix(data=c(0, Re.juv.mosq, Re.adu.mosq, Re.resjuv.mosq, Re.resadu.mosq,
                       Re.mosq.juv, 0, 0, 0, 0,
                       Re.mosq.adu, 0, 0, 0, 0,
                       Re.mosq.resjuv, 0, 0, 0, 0,
                       Re.mosq.resadu, 0, 0, 0, 0), 
                nrow=5, ncol=5)

  eigenvalues <- eigen(NGM)
  Re <- max(abs(eigenvalues$values))  
  
  print(Re)
  
  Re.juv <- sqrt(Re.mosq.juv * Re.juv.mosq)
  Re.adu <- sqrt(Re.mosq.adu * Re.adu.mosq)
  Re.bb <- sqrt(Re.mosq.juv * Re.juv.mosq + Re.mosq.adu * Re.adu.mosq) 
  Re.resjuv <- sqrt(Re.mosq.resjuv * Re.resjuv.mosq)
  Re.resadu <- sqrt(Re.mosq.resadu * Re.resadu.mosq)
  Re.res <- sqrt(Re.mosq.resjuv * Re.resjuv.mosq + Re.mosq.resadu * Re.resadu.mosq)
  
  
  Re.output <- cbind(Re.mosq.juv, Re.mosq.adu, Re.mosq.resjuv, Re.mosq.resadu, Re.juv.mosq, Re.adu.mosq, Re.resjuv.mosq, Re.resadu.mosq,
                     Re.juv, Re.adu, Re.bb, Re.resjuv, Re.resadu, Re.res, Re)
  return(Re.output)
  
}

process.R0 <- function(simOutput, average) {
  
  R0.processed <- data.frame(matrix(nrow=nrow(simOutput[[1]]), ncol=0))
  
  
  if (average=="yes") {
    
  } else if (average=="no") {
    for (ii in 1:length(simOutput)) {
      #remove last row - check what goes wrong in aggregation. solved
      # n<-dim(simOutput[[ii]])[1]
      # R0.list_i <- simOutput[[ii]][1:(n-1),]
      R0.list_i <- simOutput[[ii]]
      R0.processed <- cbind(R0.processed, R0.list_i$R0.mosq.juv, R0.list_i$R0.mosq.adu, R0.list_i$R0.mosq.resjuv, R0.list_i$R0.mosq.resadu,
                            R0.list_i$R0.juv.mosq, R0.list_i$R0.adu.mosq, R0.list_i$R0.resjuv.mosq, R0.list_i$R0.resadu.mosq, 
                            R0.list_i$R0.juv, R0.list_i$R0.adu, R0.list_i$R0.bb, R0.list_i$R0.resjuv, R0.list_i$R0.resadu, R0.list_i$R0.res,
                            R0.list_i$R0)
    }
  } 
  
  # add date column
  R0.processed <- cbind(R0.processed, R0.list_i$date)
  
  # browser()
  
  colnames(R0.processed) <- c(paste(c("R0.mosq.juv", "R0.mosq.adu", "R0.mosq.resjuv", "R0.mosq.resadu",
                                    "R0.juv.mosq", "R0.adu.mosq", "R0.resjuv.mosq", "R0.resadu.mosq",
                                    "R0.juvv", "R0.aduu", "R0.bb", "R0.resjuvv", "R0.resaduu", "R0.ress", "R0.all"), rep(1:length(simOutput), each=15), sep="_"), "date")
  

  R0.processed <- R0.processed %>%
    filter(!is.na(R0.all_1)) %>%
    rowwise(.) %>%
    mutate(mean.mosq.juv = mean(c_across(c(all_of(grep("R0.mosq.juv", names(.), fixed=TRUE)))))) %>%
    mutate(lb.mosq.juv = quantile(c_across(c(all_of(grep("R0.mosq.juv", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.mosq.juv = quantile(c_across(c(all_of(grep("R0.mosq.juv", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.mosq.adu = mean(c_across(c(all_of(grep("R0.mosq.adu", names(.), fixed=TRUE)))))) %>%
    mutate(lb.mosq.adu = quantile(c_across(c(all_of(grep("R0.mosq.adu", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.mosq.adu = quantile(c_across(c(all_of(grep("R0.mosq.adu", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.mosq.resjuv = mean(c_across(c(all_of(grep("R0.mosq.resjuv", names(.), fixed=TRUE)))))) %>%
    mutate(lb.mosq.resjuv = quantile(c_across(c(all_of(grep("R0.mosq.resjuv", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.mosq.resjuv = quantile(c_across(c(all_of(grep("R0.mosq.resjuv", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.mosq.resadu = mean(c_across(c(all_of(grep("R0.mosq.resadu", names(.), fixed=TRUE)))))) %>%
    mutate(lb.mosq.resadu = quantile(c_across(c(all_of(grep("R0.mosq.resadu", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.mosq.resadu = quantile(c_across(c(all_of(grep("R0.mosq.resadu", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%

    mutate(mean.juv.mosq = mean(c_across(c(all_of(grep("R0.juv.mosq", names(.), fixed=TRUE)))))) %>%
    mutate(lb.juv.mosq = quantile(c_across(c(all_of(grep("R0.juv.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.juv.mosq = quantile(c_across(c(all_of(grep("R0.juv.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.adu.mosq = mean(c_across(c(all_of(grep("R0.adu.mosq", names(.), fixed=TRUE)))))) %>%
    mutate(lb.adu.mosq = quantile(c_across(c(all_of(grep("R0.adu.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.adu.mosq = quantile(c_across(c(all_of(grep("R0.adu.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.resjuv.mosq = mean(c_across(c(all_of(grep("R0.resjuv.mosq", names(.), fixed=TRUE)))))) %>%
    mutate(lb.resjuv.mosq = quantile(c_across(c(all_of(grep("R0.resjuv.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.resjuv.mosq = quantile(c_across(c(all_of(grep("R0.resjuv.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.resadu.mosq = mean(c_across(c(all_of(grep("R0.resadu.mosq", names(.), fixed=TRUE)))))) %>%
    mutate(lb.resadu.mosq = quantile(c_across(c(all_of(grep("R0.resadu.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.resadu.mosq = quantile(c_across(c(all_of(grep("R0.resadu.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%

    mutate(mean.juv = mean(c_across(c(all_of(grep("R0.juvv", names(.), fixed=TRUE)))))) %>%
    mutate(lb.juv = quantile(c_across(c(all_of(grep("R0.juvv", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.juv = quantile(c_across(c(all_of(grep("R0.juvv", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.adu = mean(c_across(c(all_of(grep("R0.aduu", names(.), fixed=TRUE)))))) %>%
    mutate(lb.adu = quantile(c_across(c(all_of(grep("R0.aduu", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.adu = quantile(c_across(c(all_of(grep("R0.aduu", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%

    mutate(mean.bb = mean(c_across(c(all_of(grep("R0.bb", names(.), fixed=TRUE)))))) %>%
    mutate(lb.bb = quantile(c_across(c(all_of(grep("R0.bb", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.bb = quantile(c_across(c(all_of(grep("R0.bb", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%

    mutate(mean.resjuv = mean(c_across(c(all_of(grep("R0.resjuvv", names(.), fixed=TRUE)))))) %>%
    mutate(lb.resjuv = quantile(c_across(c(all_of(grep("R0.resjuvv", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.resjuv = quantile(c_across(c(all_of(grep("R0.resjuvv", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.resadu = mean(c_across(c(all_of(grep("R0.resaduu", names(.), fixed=TRUE)))))) %>%
    mutate(lb.resadu = quantile(c_across(c(all_of(grep("R0.resaduu", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.resadu = quantile(c_across(c(all_of(grep("R0.resaduu", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.res = mean(c_across(c(all_of(grep("R0.ress", names(.), fixed=TRUE)))))) %>%
    mutate(lb.res = quantile(c_across(c(all_of(grep("R0.ress", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.res = quantile(c_across(c(all_of(grep("R0.ress", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.all = mean(c_across(c(all_of(grep("R0.all", names(.), fixed=TRUE)))))) %>%
    mutate(lb.all = quantile(c_across(c(all_of(grep("R0.all", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.all = quantile(c_across(c(all_of(grep("R0.all", names(.), fixed=TRUE)))), c(0.025, 0.975))[2])
  
  return(R0.processed)
}

process.R0.space <- function(simOutput, average) {
  
  R0.processed <- data.frame(matrix(nrow=nrow(simOutput[[1]]), ncol=0))
  
  
  if (average=="yes") {
    
  } else if (average=="no") {
    for (ii in 1:length(simOutput)) {
      R0.list_i <- simOutput[[ii]]
      R0.processed <- cbind(R0.processed, R0.list_i$R0)
    }
  } 
  
  # add node and year column
  R0.processed <- cbind(R0.processed, R0.list_i$node, R0.list_i$year)
  
  # browser()
  colnames(R0.processed) <- c(paste(c("R0.all"), rep(1:length(simOutput), each=1), sep="_"), "node", "year")
  
  R0.processed <- R0.processed %>%
    filter(!is.na(R0.all_1)) %>%
    rowwise(.) %>%
    mutate(mean.all = mean(c_across(c(all_of(grep("R0.all", names(.), fixed=TRUE)))))) %>%
    mutate(lb.all = quantile(c_across(c(all_of(grep("R0.all", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.all = quantile(c_across(c(all_of(grep("R0.all", names(.), fixed=TRUE)))), c(0.025, 0.975))[2])
  
  return(R0.processed)
}

process.Re <- function(simOutput, average) {
  
  Re.processed <- data.frame(matrix(nrow=nrow(simOutput[[1]]), ncol=0))
  
  
  if (average=="yes") {
 
  } else if (average=="no") {
    for (ii in 1:length(simOutput)) {

      Re.list_i <- simOutput[[ii]]
      Re.processed <- cbind(Re.processed, Re.list_i$Re.mosq.juv, Re.list_i$Re.mosq.adu, Re.list_i$Re.mosq.resjuv, Re.list_i$Re.mosq.resadu,
                            Re.list_i$Re.juv.mosq, Re.list_i$Re.adu.mosq, Re.list_i$Re.resjuv.mosq, Re.list_i$Re.resadu.mosq, 
                            Re.list_i$Re.juv, Re.list_i$Re.adu, Re.list_i$Re.bb, Re.list_i$Re.resjuv, Re.list_i$Re.resadu, Re.list_i$Re.res,
                            Re.list_i$Re)    }
  }
  
  # add date column
  Re.processed <- cbind(Re.processed, Re.list_i$date)
  
  # browser()
  
  colnames(Re.processed) <- c(paste(c("Re.mosq.juv", "Re.mosq.adu", "Re.mosq.resjuv", "Re.mosq.resadu",
                                      "Re.juv.mosq", "Re.adu.mosq", "Re.resjuv.mosq", "Re.resadu.mosq",
                                      "Re.juvv", "Re.aduu", "Re.bb", "Re.resjuvv", "Re.resaduu", "Re.ress", "Re.all"), rep(1:length(simOutput), each=15), sep="_"), "date")
  
  Re.processed <- Re.processed %>%
    filter(!is.na(Re.all_1)) %>%
    rowwise(.) %>%
    mutate(mean.mosq.juv = mean(c_across(c(all_of(grep("Re.mosq.juv", names(.), fixed=TRUE)))))) %>%
    mutate(lb.mosq.juv = quantile(c_across(c(all_of(grep("Re.mosq.juv", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.mosq.juv = quantile(c_across(c(all_of(grep("Re.mosq.juv", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.mosq.adu = mean(c_across(c(all_of(grep("Re.mosq.adu", names(.), fixed=TRUE)))))) %>%
    mutate(lb.mosq.adu = quantile(c_across(c(all_of(grep("Re.mosq.adu", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.mosq.adu = quantile(c_across(c(all_of(grep("Re.mosq.adu", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.mosq.resjuv = mean(c_across(c(all_of(grep("Re.mosq.resjuv", names(.), fixed=TRUE)))))) %>%
    mutate(lb.mosq.resjuv = quantile(c_across(c(all_of(grep("Re.mosq.resjuv", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.mosq.resjuv = quantile(c_across(c(all_of(grep("Re.mosq.resjuv", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.mosq.resadu = mean(c_across(c(all_of(grep("Re.mosq.resadu", names(.), fixed=TRUE)))))) %>%
    mutate(lb.mosq.resadu = quantile(c_across(c(all_of(grep("Re.mosq.resadu", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.mosq.resadu = quantile(c_across(c(all_of(grep("Re.mosq.resadu", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.juv.mosq = mean(c_across(c(all_of(grep("Re.juv.mosq", names(.), fixed=TRUE)))))) %>%
    mutate(lb.juv.mosq = quantile(c_across(c(all_of(grep("Re.juv.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.juv.mosq = quantile(c_across(c(all_of(grep("Re.juv.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.adu.mosq = mean(c_across(c(all_of(grep("Re.adu.mosq", names(.), fixed=TRUE)))))) %>%
    mutate(lb.adu.mosq = quantile(c_across(c(all_of(grep("Re.adu.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.adu.mosq = quantile(c_across(c(all_of(grep("Re.adu.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.resjuv.mosq = mean(c_across(c(all_of(grep("Re.resjuv.mosq", names(.), fixed=TRUE)))))) %>%
    mutate(lb.resjuv.mosq = quantile(c_across(c(all_of(grep("Re.resjuv.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.resjuv.mosq = quantile(c_across(c(all_of(grep("Re.resjuv.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.resadu.mosq = mean(c_across(c(all_of(grep("Re.resadu.mosq", names(.), fixed=TRUE)))))) %>%
    mutate(lb.resadu.mosq = quantile(c_across(c(all_of(grep("Re.resadu.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.resadu.mosq = quantile(c_across(c(all_of(grep("Re.resadu.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.juv = mean(c_across(c(all_of(grep("Re.juvv", names(.), fixed=TRUE)))))) %>%
    mutate(lb.juv = quantile(c_across(c(all_of(grep("Re.juvv", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.juv = quantile(c_across(c(all_of(grep("Re.juvv", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.adu = mean(c_across(c(all_of(grep("Re.aduu", names(.), fixed=TRUE)))))) %>%
    mutate(lb.adu = quantile(c_across(c(all_of(grep("Re.aduu", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.adu = quantile(c_across(c(all_of(grep("Re.aduu", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.bb = mean(c_across(c(all_of(grep("Re.bb", names(.), fixed=TRUE)))))) %>%
    mutate(lb.bb = quantile(c_across(c(all_of(grep("Re.bb", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.bb = quantile(c_across(c(all_of(grep("Re.bb", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.resjuv = mean(c_across(c(all_of(grep("Re.resjuvv", names(.), fixed=TRUE)))))) %>%
    mutate(lb.resjuv = quantile(c_across(c(all_of(grep("Re.resjuvv", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.resjuv = quantile(c_across(c(all_of(grep("Re.resjuvv", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.resadu = mean(c_across(c(all_of(grep("Re.resaduu", names(.), fixed=TRUE)))))) %>%
    mutate(lb.resadu = quantile(c_across(c(all_of(grep("Re.resaduu", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.resadu = quantile(c_across(c(all_of(grep("Re.resaduu", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.res = mean(c_across(c(all_of(grep("Re.ress", names(.), fixed=TRUE)))))) %>%
    mutate(lb.res = quantile(c_across(c(all_of(grep("Re.ress", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.res = quantile(c_across(c(all_of(grep("Re.ress", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.all = mean(c_across(c(all_of(grep("Re.all", names(.), fixed=TRUE)))))) %>%
    mutate(lb.all = quantile(c_across(c(all_of(grep("Re.all", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.all = quantile(c_across(c(all_of(grep("Re.all", names(.), fixed=TRUE)))), c(0.025, 0.975))[2])
  
  return(Re.processed)
}

process.Re.space <- function(simOutput, average) {
  
  Re.processed <- data.frame(matrix(nrow=nrow(simOutput[[1]]), ncol=0))
  
  if (average=="yes") {

  } else if (average=="no") {
    for (ii in 1:length(simOutput)) {
      Re.list_i <- simOutput[[ii]]
      Re.processed <- cbind(Re.processed, Re.list_i$Re)    }
  }
  
  # add date column
  Re.processed <- cbind(Re.processed, Re.list_i$date, Re.list_i$year)
  # browser()
  colnames(Re.processed) <- c(paste(c("Re.all"), rep(1:length(simOutput), each=1), sep="_"), "node", "year")
  
  Re.processed <- Re.processed %>%
    filter(!is.na(Re.all_1)) %>%
    rowwise(.) %>%
    mutate(mean.all = mean(c_across(c(all_of(grep("Re.all", names(.), fixed=TRUE)))))) %>%
    mutate(lb.all = quantile(c_across(c(all_of(grep("Re.all", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.all = quantile(c_across(c(all_of(grep("Re.all", names(.), fixed=TRUE)))), c(0.025, 0.975))[2])
  
  return(Re.processed)
}

# Prevalence --------------------------------------------------------------

calculate.prev <- function(output) {
  
  # browser()
  Njuv <- output$Shj + output$Ehj + output$Ihj + output$Rhj
  Nadu <- output$Sha + output$Eha + output$Iha + output$Rha
  Nmosq <- output$Sm + output$Em + output$Im
  
  prev.Juv <- (output$Ehj + output$Ihj) / Njuv
  prev.Adu <- (output$Eha + output$Iha) / Nadu
  prev.both <- (output$Ehj + output$Ihj + output$Eha + output$Iha) / (Njuv + Nadu)
  prev.mosq <- (output$Em + output$Im) / Nmosq
  
  prev <- cbind(prev.Juv, prev.Adu, prev.both, prev.mosq)
  
  return(prev)
}

process.prev <- function(simOutput, average) {
  
  prev.processed <- data.frame(matrix(nrow=nrow(simOutput[[1]]), ncol=0))
  
  
  if (average=="yes") {
    
  } else if (average=="no") {
    for (ii in 1:length(simOutput)) {
      prev.list_i <- simOutput[[ii]]
      prev.processed <- cbind(prev.processed, prev.list_i$prev.Juv, prev.list_i$prev.Adu, prev.list_i$prev.both, prev.list_i$prev.mosq)
    }
  }
  
  # add date column
  prev.processed <- cbind(prev.processed, prev.list_i$date)
  
  # browser()
  colnames(prev.processed) <- c(paste(c("prev.Juv", "prev.Adu", "prev.both", "prev.mosq"), 
                                      rep(1:length(simOutput), each = 4), sep="_"), "date")
  
  prev.processed <- prev.processed %>%
    filter(!is.na(prev.both_1)) %>%
    rowwise(.) %>%
    mutate(mean.Juv = mean(c_across(c(all_of(grep("prev.Juv", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.Juv = quantile(c_across(c(all_of(grep("prev.Juv", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.Juv = quantile(c_across(c(all_of(grep("prev.Juv", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) %>%
    
    mutate(mean.Adu = mean(c_across(c(all_of(grep("prev.Adu", names(.), fixed=TRUE)))))) %>%
    mutate(lb.Adu = quantile(c_across(c(all_of(grep("prev.Adu", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.Adu = quantile(c_across(c(all_of(grep("prev.Adu", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.both = mean(c_across(c(all_of(grep("prev.both", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.both = quantile(c_across(c(all_of(grep("prev.both", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.both = quantile(c_across(c(all_of(grep("prev.both", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) %>%

    mutate(mean.mosq = mean(c_across(c(all_of(grep("prev.mosq", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.mosq = quantile(c_across(c(all_of(grep("prev.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.mosq = quantile(c_across(c(all_of(grep("prev.mosq", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) 
  
  return(prev.processed)
}

# Sero prevalence --------------------------------------------------------------

calculate.sero <- function(output) {
  
  # browser()
  Njuv <- output$Shj + output$Ehj + output$Ihj + output$Rhj
  Nadu <- output$Sha + output$Eha + output$Iha + output$Rha
  Nresjuv <- output$Srj + output$Erj + output$Irj + output$Rrj
  Nresadu <- output$Sra + output$Era + output$Ira + output$Rra
  
  sero.Juv <- output$Rhj / Njuv
  sero.Adu <- output$Rha / Nadu
  sero.bothbb <- (output$Rhj + output$Rha) / (Njuv + Nadu)
  sero.resjuv <- output$Rrj / Nresjuv
  sero.resadu <- output$Rra / Nresadu
  sero.bothres <- (output$Rrj + output$Rra) / (Nresjuv + Nresadu)
  
  sero <- cbind(sero.Juv, sero.Adu, sero.bothbb, sero.resjuv, sero.resadu, sero.bothres)
  
  return(sero)
}

process.sero <- function(simOutput, average) {
  
  sero.processed <- data.frame(matrix(nrow=nrow(simOutput[[1]]), ncol=0))
  
  
  if (average=="yes") {

  } else if (average=="no") {
    for (ii in 1:length(simOutput)) {
      sero.list_i <- simOutput[[ii]]
      sero.processed <- cbind(sero.processed, sero.list_i$sero.Juv, sero.list_i$sero.Adu, sero.list_i$sero.bothbb, 
                              sero.list_i$sero.resjuv, sero.list_i$sero.resadu, sero.list_i$sero.bothres)
    }
  }
  
  # add date column
  sero.processed <- cbind(sero.processed, sero.list_i$date)
  
  # browser()
  colnames(sero.processed) <- c(paste(c("sero.Juv", "sero.Adu", "sero.bothbb", "sero.resjuv", "sero.resadu", "sero.bothres"), rep(1:length(simOutput), each = 6), sep="_"), "date")
  
  sero.processed <- sero.processed %>%
    filter(!is.na(sero.bothbb_1)) %>%
    rowwise(.) %>%
    mutate(mean.Juv = mean(c_across(c(all_of(grep("sero.Juv", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.Juv = quantile(c_across(c(all_of(grep("sero.Juv", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.Juv = quantile(c_across(c(all_of(grep("sero.Juv", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) %>%
    
    mutate(mean.Adu = mean(c_across(c(all_of(grep("sero.Adu", names(.), fixed=TRUE)))))) %>%
    mutate(lb.Adu = quantile(c_across(c(all_of(grep("sero.Adu", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.Adu = quantile(c_across(c(all_of(grep("sero.Adu", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.bothbb = mean(c_across(c(all_of(grep("sero.bothbb", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.bothbb = quantile(c_across(c(all_of(grep("sero.bothbb", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.bothbb = quantile(c_across(c(all_of(grep("sero.bothbb", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) %>%

    mutate(mean.resjuv = mean(c_across(c(all_of(grep("sero.resjuv", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.resjuv = quantile(c_across(c(all_of(grep("sero.resjuv", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.resjuv = quantile(c_across(c(all_of(grep("sero.resjuv", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) %>%
    
    mutate(mean.resadu = mean(c_across(c(all_of(grep("sero.resadu", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.resadu = quantile(c_across(c(all_of(grep("sero.resadu", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.resadu = quantile(c_across(c(all_of(grep("sero.resadu", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) %>%
    
    mutate(mean.bothres = mean(c_across(c(all_of(grep("sero.bothres", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.bothres = quantile(c_across(c(all_of(grep("sero.bothres", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.bothres = quantile(c_across(c(all_of(grep("sero.bothres", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) 
  
  return(sero.processed)
}

# Dead prevalence --------------------------------------------------------------

calculate.newdeaths <- function(output) {
  # browser()
  
  output <- output %>%
    arrange(date) %>%
    mutate(Ddhj_new = Ddhj - lag(Ddhj, n=1)) %>%
    mutate(DIhj_new = DIhj - lag(DIhj, n=1)) %>%
    mutate(Ddha_new = Ddha - lag(Ddha, n=1)) %>%
    mutate(DIha_new = DIha - lag(DIha, n=1))
  
  return(output)
  
}

calculate.dead <- function(output) {
    # browser()

  dead.Juv <- output$DIhj_new / (output$DIhj_new + output$Ddhj_new)
  dead.Adu <- output$DIha_new / (output$DIha_new + output$Ddha_new)
  dead.both <- (output$DIhj_new + output$DIha_new) / (output$DIhj_new + output$Ddhj_new + output$DIha_new + output$Ddha_new)

  dead <- cbind(dead.Juv, dead.Adu, dead.both)
  
  return(dead)
}

process.dead <- function(simOutput, average) {
  
  dead.processed <- data.frame(matrix(nrow=nrow(simOutput[[1]]), ncol=0))
  
  
  if (average=="yes") {

      } else if (average=="no") {
    for (ii in 1:length(simOutput)) {

      dead.list_i <- simOutput[[ii]]
      dead.processed <- cbind(dead.processed, dead.list_i$dead.Juv, dead.list_i$dead.Adu, dead.list_i$dead.both)
    }
  }
  
  # add date column
  dead.processed <- cbind(dead.processed, dead.list_i$date)
  
  # browser()
  colnames(dead.processed) <- c(paste(c("dead.Juv", "dead.Adu", "dead.both"), rep(1:length(simOutput), each = 3), sep="_"), "date")
  
  dead.processed <- dead.processed %>%
    filter(!is.na(dead.both_1)) %>%
    rowwise(.) %>%
    mutate(mean.Juv = mean(c_across(c(all_of(grep("dead.Juv", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.Juv = quantile(c_across(c(all_of(grep("dead.Juv", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.Juv = quantile(c_across(c(all_of(grep("dead.Juv", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) %>%
    
    mutate(mean.Adu = mean(c_across(c(all_of(grep("dead.Adu", names(.), fixed=TRUE)))))) %>%
    mutate(lb.Adu = quantile(c_across(c(all_of(grep("dead.Adu", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.Adu = quantile(c_across(c(all_of(grep("dead.Adu", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
    
    mutate(mean.both = mean(c_across(c(all_of(grep("dead.both", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.both = quantile(c_across(c(all_of(grep("dead.both", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.both = quantile(c_across(c(all_of(grep("dead.both", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) 
  
  return(dead.processed)
}


# Biting preference -------------------------------------------------------

calculate.biting.bb <- function(output, posterior) {
  
  bitingPref.parameter <- posterior$bitingPref
  Njuv <- output$Shj + output$Ehj + output$Ihj + output$Rhj
  Nadu <- output$Sha + output$Eha + output$Iha + output$Rha

  Nbb <- Njuv + Nadu
  Nres <- output$Srj + output$Erj + output$Irj + output$Rrj + output$Sra + output$Era + output$Ira + output$Rra

  
  prop.bites.blackbird.vs.res <- Nbb / (Nbb + bitingPref.parameter*Nres)
  prop.bites.blacbird.vs.total <- prop.bites.blackbird.vs.res * output$prop.bites
  res.vsbb <- bitingPref.parameter*Nres / Nbb
  
  bitingPref <- cbind(Nbb, Nres, prop.bites.blackbird.vs.res, prop.bites.blacbird.vs.total, res.vsbb)
  
  return(bitingPref)
  
}

process.biting <- function(simOutput) {
  
  # browser()
  biting.processed <- data.frame(matrix(nrow=nrow(simOutput[[1]]), ncol=0))
  
    for (ii in 1:length(simOutput)) {
      biting.list_i <- simOutput[[ii]]
      biting.processed <- cbind(biting.processed, biting.list_i$Nbb, biting.list_i$Nres, biting.list_i$prop.bites,
                                biting.list_i$prop.bites.bb.vs.res, 
                                biting.list_i$prop.bites.bb.vs.total,
                                biting.list_i$res.vsbb)
    }
  
  # add date column
  biting.processed <- cbind(biting.processed, biting.list_i$date)
  
  # browser()
  colnames(biting.processed) <- c(paste(c("Nbb", "Nres", "propBites", "prop.bites.blackbird.vs.res", "prop.bites.blackbird.vs.total", "res.vsbb"), 
                                        rep(1:length(simOutput), each = 6), sep="_"), "date")
  
  biting.processed <- biting.processed %>%
    filter(!is.na(prop.bites.blackbird.vs.res_1)) %>%
    rowwise(.) %>%
    mutate(mean.Nbb = mean(c_across(c(all_of(grep("Nbb", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.Nbb = quantile(c_across(c(all_of(grep("Nbb", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.Nbb = quantile(c_across(c(all_of(grep("Nbb", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) %>%
    
    mutate(mean.Nres = mean(c_across(c(all_of(grep("Nres", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.Nres = quantile(c_across(c(all_of(grep("Nres", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.Nres = quantile(c_across(c(all_of(grep("Nres", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) %>%
    
    mutate(mean.propbites = mean(c_across(c(all_of(grep("propBites", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.propbites = quantile(c_across(c(all_of(grep("propBites", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.propbites = quantile(c_across(c(all_of(grep("propBites", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) %>%
    
    mutate(mean.vsres = mean(c_across(c(all_of(grep("prop.bites.blackbird.vs.res", names(.), fixed=TRUE)))), na.rm=T)) %>%
    mutate(lb.vsres = quantile(c_across(c(all_of(grep("prop.bites.blackbird.vs.res", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[1]) %>%
    mutate(ub.vsres = quantile(c_across(c(all_of(grep("prop.bites.blackbird.vs.res", names(.), fixed=TRUE)))), c(0.025, 0.975), na.rm=T)[2]) %>%
    
    mutate(mean.vstot = mean(c_across(c(all_of(grep("prop.bites.blackbird.vs.total", names(.), fixed=TRUE)))))) %>%
    mutate(lb.vstot = quantile(c_across(c(all_of(grep("prop.bites.blackbird.vs.total", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.vstot = quantile(c_across(c(all_of(grep("prop.bites.blackbird.vs.total", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) %>%
  
    mutate(mean.resvsbb = mean(c_across(c(all_of(grep("res.vsbb", names(.), fixed=TRUE)))))) %>%
    mutate(lb.resvsbb = quantile(c_across(c(all_of(grep("res.vsbb", names(.), fixed=TRUE)))), c(0.025, 0.975))[1]) %>%
    mutate(ub.resvsbb = quantile(c_across(c(all_of(grep("res.vsbb", names(.), fixed=TRUE)))), c(0.025, 0.975))[2]) 
    
  return(biting.processed)
}

