#________________________________________________________________________________

## Description:
# This files contains the parameters for the analysis/paper called: 

#________________________________________________________________________________

load("../Data/model_input/temperature5kmNoWater.RData")

load("../Data/model_input/dispersal_daily.bb.RData")
load("../Data/model_input/dispersal_daily.sc20.RData")


# Settings ----------------------------------------------------------------

n <- nrow(temperature)

# Compartments & transitions ----------------------------------------------

compartments <- c("Sm", "Em", "Im", 
                  "Shj", "Ehj", "Ihj", "Rhj", "Duhj", "Ddhj", "DIhj", 
                  "Sha", "Eha", "Iha", "Rha", "Duha", "Ddha", "DIha",
                  "Srj", "Erj", "Irj", "Rrj", "Drj",
                  "Sra", "Era", "Ira", "Rra", "Dra")

transitions <- c("Nhj <- Shj + Ehj + Ihj + Rhj",
                 "Nha <- Sha + Eha + Iha + Rha",
                 "Nrj <- Srj + Erj + Irj + Rrj",
                 "Nra <- Sra + Era + Ira + Rra",
  
                 "Sm -> biting * prop.bites * ((Nhj+0.001) / ((Nrj+Nra)* pref + (Nhj+Nha+0.001))) * transmissionProbHM * Sm *  I_coupling_juv -> Em",
                 "Sm -> biting * prop.bites * ((Nha+0.001) / ((Nrj+Nra)* pref + (Nhj+Nha+0.001))) * transmissionProbHM * Sm  * I_coupling_adu -> Em",
                 "Sm -> biting * prop.bites * ((Nrj+0.001)*pref / ((Nrj+Nra)* pref + (Nhj+Nha+0.001))) * transmissionProbHM * Sm  * I_coupling_resjuv -> Em",
                 "Sm -> biting * prop.bites * ((Nra+0.001)*pref / ((Nrj+Nra)* pref + (Nhj+Nha+0.001))) * transmissionProbHM * Sm  * I_coupling_resadu -> Em",
                 "Sm -> reemergence * Sm -> Em",
                 "Em -> incubationRateM * Em -> Im",
                 "Sm -> diapause * Sm -> @",
                 "Em -> diapause * Em -> @",
                 "Im -> diapause * Im -> @",

                 "@ -> birthRateH * Nha -> Shj",
                 "Shj -> biting * prop.bites * ((Nhj+0.001) / ((Nrj+Nra)* pref + (Nhj+Nha+0.001))) * transmissionProbMH * Shj  * S_coupling_juv -> Ehj",
                 "Ehj -> incubationRateH * Ehj -> Ihj",
                 "Ihj -> recoveryRateH * Ihj -> Rhj",
                 "Ihj -> disIndMortalityH * Ihj -> DIhj",
                 "Rhj -> waning * Rhj -> Shj",
                 "Shj -> mortalityRateHJuv * Shj * (1-detection) -> Duhj",
                 "Shj -> mortalityRateHJuv * Shj * detection -> Ddhj",
                 "Ehj -> mortalityRateHJuv * Ehj -> DIhj",
                 "Ihj -> mortalityRateHJuv * Ihj -> DIhj",
                 "Rhj -> mortalityRateHJuv * Rhj * (1-detection) -> Duhj",
                 "Rhj -> mortalityRateHJuv * Rhj * detection -> Ddhj",

                 "Sha -> biting * prop.bites * ((Nha+0.001) / ((Nrj+Nra)* pref + (Nhj+Nha+0.001))) * transmissionProbMH * Sha *  S_coupling_adu -> Eha",
                 "Eha -> incubationRateH * Eha -> Iha",
                 "Iha -> recoveryRateH * Iha -> Rha",
                 "Iha -> disIndMortalityH * Iha -> DIha",
                 "Rha -> waning * Rha -> Sha",
                 "Sha -> mortalityRateHAdu * Sha * (1-detection) -> Duha",
                 "Sha -> mortalityRateHAdu * Sha * detection -> Ddha",
                 "Eha -> mortalityRateHAdu * Eha -> DIha",
                 "Iha -> mortalityRateHAdu * Iha -> DIha",
                 "Rha -> mortalityRateHAdu * Rha * (1-detection) -> Duha",
                 "Rha -> mortalityRateHAdu * Rha * detection -> Ddha",

                 "@ -> birthRateR * Nra -> Srj",
                 "Srj -> biting * prop.bites * ((Nrj+0.001)*pref / ((Nrj+Nra)* pref + (Nhj+Nha+0.001))) * transmissionProbMR * Srj * S_coupling_resjuv -> Erj",
                 "Erj -> incubationRateH * Erj -> Irj",
                 "Irj -> recoveryRateH * Irj -> Rrj",
                 "Irj -> disIndMortalityR * Irj -> Drj",
                 "Rrj -> waning * Rrj -> Srj",
                 "Srj -> mortalityRateHJuv * Srj -> Drj",
                 "Erj -> mortalityRateHJuv * Erj -> Drj",
                 "Irj -> mortalityRateHJuv * Irj -> Drj",
                 "Rrj -> mortalityRateHJuv * Rrj -> Drj",
                 
                 "Sra -> biting * prop.bites * ((Nrj+0.001)*pref / ((Nrj+Nra)* pref + (Nhj+Nha+0.001))) * transmissionProbMR * Sra * S_coupling_resadu -> Era",
                 "Era -> incubationRateH * Era -> Ira",
                 "Ira -> recoveryRateH * Ira -> Rra",
                 "Ira -> disIndMortalityR * Ira -> Dra",
                 "Rra -> waning * Rra -> Sra",
                 "Sra -> (1/lifespanR) * Sra -> Dra",
                 "Era -> (1/lifespanR) * Era -> Dra",
                 "Ira -> (1/lifespanR) * Ira -> Dra",
                 "Rra -> (1/lifespanR) * Rra -> Dra")



# Parameters --------------------------------------------------------------

# parameters outside of ldata, to be used in u0
scalingParameter <- 0.25 

FOI.n <- rep(0.13, n/3)
FOI.m <- rep(0.13, n/3)
FOI.s <- rep(0.13, n/3)
FOI <- c(FOI.n, FOI.m, FOI.s)

ldata <- data.frame(
  bitingRate1 = rep(0.344, n),        # biting rate: 0.344/(1+1.231*exp(-0.184*(temp-20)))
  bitingRate2 = rep(1, n),
  bitingRate3 = rep(1.231, n),
  bitingRate4 = rep(0.184, n),
  bitingRate5 = rep(20, n),
  transmissionProbHM = rep(0.74, n),   
  transmissionProbMH = rep(0.88, n),  
  transmissionProbMR = rep(0.88, n),  
  detection = rep(1, n),              # not used in current version, but allows to specify a detection probability for dead birds based on infection status
  birthRateH = rep(0.0, n),           # birth rate bird included through 'events' (see functions_run_model_x.R)
  birthRateR = rep(0, n),             # birth rate bird included through 'events' (see functions_run_model_x.R)
  mortalityRateHAdu = rep(0.0011, n), 
  mortalityRateHJuv = rep(0.0059, n), 
  lifespanR = rep(909, n),            
  incubationRateH = rep(0.67, n),      
  recoveryRateH = rep(0.25, n),        
  recoveryRateR = rep(0.25, n),        
  disIndMortalityH = rep(0.75, n),     
  disIndMortalityR = rep(0, n),     
  incubationRateM1 = rep(7.38*10^-5, n),  # extrinsic incubation rate // For T>13◦C: 7.38*10^-5*T*(T-11.4)*(45.2-T)^0.5) else 0. Shocket paper
  incubationRateM2 = rep(11.4, n), 
  incubationRateM3 = rep(45.2,n), 
  mortalityRateM1 = rep(4.86,n),      # death rate mosquitoes // For T<32: 1/(-4.86*T + 169.8), else 0.17 Shocket et al
  mortalityRateM2 = rep(169.8,n),
  mortalityRateMHighTemp = rep(0.17,n),
  temperatureCutoff1 = rep(13,n),     # used in functional form of incubationRateM
  temperatureCutoff2 = rep(32,n),     # used in functional form of mortalityRateM
  vertical_transm = rep(0.19),         
  diapauseRate = rep(0.05,n),         
  waning = rep(0.0014,n),              
  pref = rep(7.6,n)  
) 

# Creation ldata ----------------------------------------------------------
# combined temperature, parameters and dispersal probabilities

# take out x, y, layer columns 
temperature <- temperature[,-c(1:3)] 

# bind temperature and dispersal data to ldata
# dispersal matrix is different for each dispersal scenario, create one ldata file per scenario
# first dispersal is for main host (always bb), second set is for reservoir host
ldata_dispersal.bb = cbind(temperature, ldata, dispersal_daily.bb[[1]], dispersal_daily.bb[[2]], dispersal_daily.bb[[1]], dispersal_daily.bb[[2]])
ldata_dispersal.sc20 = cbind(temperature, ldata, dispersal_daily.bb[[1]], dispersal_daily.bb[[2]], dispersal_daily.sc20[[1]], dispersal_daily.sc20[[2]])
ldata_dispersal.sens = cbind(temperature, ldata, dispersal_daily.bb[[1]], dispersal_daily.bb[[1]], dispersal_daily.bb[[1]], dispersal_daily.bb[[1]])
rm(ldata, temperature, dispersal_daily.bb, dispersal_daily.sc20) 


# duplicates in colnames for version: ldata_dispersal.bb. This happens because the same dispersal dataframe is used twice -> alter names
oldnames <- colnames(ldata_dispersal.bb)
first <- oldnames[1:(ncol(ldata_dispersal.bb)-2*n)] # select names of all columns up to reservoir-related movement columns (these names can stay)
second <- oldnames[(ncol(ldata_dispersal.bb)-2*n+1):ncol(ldata_dispersal.bb)] # select names of reservoir-related movement columns
second <- gsub("bb", "res", second) # change part of the column names for reservoir
colnames(ldata_dispersal.bb) <- c(first, second)
rm(first, second, oldnames)

# duplicates in colnames for version: ldata_dispersal.sens. This happens because the same dispersal dataframe is used 4 times -> alter names
oldnames <- colnames(ldata_dispersal.sens)
bb.breeding <- oldnames[1:(ncol(ldata_dispersal.sens)-3*n)] # select names of all columns up to blackbird summer-related movement columns (these names can stay)
bb.summer <- oldnames[(ncol(ldata_dispersal.sens)-3*n+1):(ncol(ldata_dispersal.sens)-2*n)] # select names of blackbird summer-related movement columns
res.breeding <- oldnames[(ncol(ldata_dispersal.sens)-2*n+1):(ncol(ldata_dispersal.sens)-1*n)] # select names of reservoir breeding season movements columns
res.summer <- oldnames[(ncol(ldata_dispersal.sens)-n+1):(ncol(ldata_dispersal.sens))] # select names of reservoir summer movements columns
bb.summer <- gsub("br", "sum", bb.summer)
res.breeding <- gsub("bb", "res", res.breeding)
res.summer <- gsub("bb_br", "res_sum", res.summer)
colnames(ldata_dispersal.sens) <- c(bb.breeding, bb.summer, res.breeding, res.summer)
rm(bb.breeding, bb.summer, res.breeding, res.summer)

# v0 -----------------------------------------------------------------------
# Starting values for temperature-dependent parameters

v0=data.frame(biting=rep(0.2,n), incubationRateM=rep(0.01,n), mumos=rep(0.01,n), 
              diapause=rep(0,n), 
              I_coupling_juv=rep(0,n), S_coupling_juv=rep(0,n),
              I_coupling_adu=rep(0,n), S_coupling_adu=rep(0,n),
              prev_aug=rep(0,n), reemergence=rep(0,n), birthRate=rep(0,n), totalBirds=rep(0,n),
              maxBirds=rep(0,n), prop.bites=rep(0,n),
              I_coupling_resjuv=rep(0,n), S_coupling_resjuv=rep(0,n),
              I_coupling_resadu=rep(0,n), S_coupling_resadu=rep(0,n))



# E matrix (for events) ---------------------------------------------------

# E matrix of size #compartments * #event types

E <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0),
            nrow = 27, ncol = 10,
            dimnames = list(compartments, c("Mosq emergence", "Mosq death", "Mosq diapause E", "Mosq diapause I", 
                                            "Infectious introductions", "Juv blackbirds", "Adu blackbirds", "Blackbird birth", 
                                            "Reservoir birth", "Juv reservoir")))

N <- matrix(c(0, 0, 0, 7, 7, 7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5, 0, 0, 0, 0, 0, 0), 
            nrow = 27, ncol = 4,
            dimnames = list(compartments, c("Ageing blackbird", "Dispersal", "Seeding", "Ageing reservoir")))


# Post time step function ---------------------------------------------------

pts_fun <- "
      double temp = ldata[(int)t-1];
      double bitingRate1 = ldata[2557];
      double bitingRate2 = ldata[2558];
      double bitingRate3 = ldata[2559];
      double bitingRate4 = ldata[2560];
      double bitingRate5 = ldata[2561];
      double birthRateR = ldata[2567];
      double incubationRateM1 = ldata[2576];
      double incubationRateM2 = ldata[2577];
      double incubationRateM3 = ldata[2578];
      double mortalityRateM1 = ldata[2579];
      double mortalityRateM2 = ldata[2580];
      double mortalityRateMHighTemp = ldata[2581];
      double temperatureCutoff1 = ldata[2582];
      double temperatureCutoff2 = ldata[2583];
      double vert_transm = ldata[2584];
      double diapauseRate = ldata[2585];

      /* Set biting rate at 0 during winter. To avoid any remaining mosquitoes to get infected */

      if (((t>289) & (t<457)) | ((t>654) & (t<822)) | ((t>1019) & (t<1187)) | ((t>1384) & (t<1552)) | 
          ((t>1750) & (t<1918)) | ((t>2115) & (t<2283)) | ((t>2480) & (t<2648)) ) {
        v_new[0] = 0;
      } else {
        v_new[0] = bitingRate1/(bitingRate2+bitingRate3*exp(-bitingRate4*(temp-bitingRate5)));
      }

      if (temp>temperatureCutoff1) {
          v_new[1] = incubationRateM1*temp*(temp-incubationRateM2)*sqrt(incubationRateM3-temp);
      } else if (temp<=temperatureCutoff1) {
          v_new[1] = 0;
      }

      if (temp<=temperatureCutoff2) {
          v_new[2] = 1/(-mortalityRateM1*temp+mortalityRateM2);
      } else if (temp>temperatureCutoff2) {
          v_new[2] = mortalityRateMHighTemp;
      }

      /* Mosquito mortality not used from Sept 15 onwards as this gets all absorbed by diapause rate. (Mortality parameter only used in NGM - mortality included in events) */
      if (((t>259) & (t<457)) | ((t>624) & (t<822)) | ((t>989) & (t<1187)) | ((t>1354) & (t<1552)) | 
          ((t>1720) & (t<1918)) | ((t>2085) & (t<2283)) | ((t>2450) & (t<2648)) ) {
          v_new[2] = 0;
      } else {
          v_new[2] = v_new[2];
      }

      if (((t>244) & (t<290)) | ((t>609) & (t<655)) | ((t>974) & (t<1020)) | ((t>1339) & (t<1385)) | 
          ((t>1705) & (t<1751)) | ((t>2070) & (t<2116)) | ((t>2435) & (t<2481)) ) {
          v_new[3] = diapauseRate;
      } else {
          v_new[3] = 0;
      }

      
      /* Determine prevalence in mosquitoes on Sept 15.
       * This, times vertical transmission, is the prevalence in mosquitoes when they leave diapause. */
      if ((t==259) | (t==624) | (t==989) | (t==1354) | (t==1720) | (t==2085)) {
        double mosq_S = u[0];
        double mosq_E = u[1];
        double mosq_I = u[2];
        double mosq_prev = mosq_I/(mosq_S + mosq_E + mosq_I);
        v_new[8] = mosq_prev;
      } else {
        v_new[8] = v[8];
      }

      if ((t==458) | (t==823) | (t==1188) | (t==1554) | (t==1919) | (t==2284) ) {
        v_new[9] = v_new[8] * vert_transm;
      } else {
        v_new[9] = 0;
      }

      /* Seasonal birth rate reservoir birds. Not used anymore, now through events. */
      if (((t>121) & (t<168)) | ((t>486) & (t<533)) | ((t>851) & (t<898)) | ((t>1216) & (t<1263)) | 
          ((t>1582) & (t<1629)) | ((t>1947) & (t<1994)) | ((t>2312) & (t<2359)) ) {
        v_new[10] = 0;
      } else {
        v_new[10] = 0;
      }

      /* Set prop.bites at each timestep and node */
      double bird_tot = u[3] + u[4] + u[5] + u[6] + u[10] + u[11] + u[12] + u[13] + u[17] + u[18] + u[19] + u[20] + u[22] + u[23] + u[24] + u[25];     
      v_new[11] = bird_tot;
      
      /* Number of compartments in the model. */
      const int Nc = 27;
      
      /* Number of nodes in the model. Nnodes=1793 with water included */
      const int Nnodes = 1398;
      
      /* Obtain matrix with current state variables at each node */
      const int *u_cur = &u[-Nc*node];

      /* Find maximum number of birds at each time step by looping over all nodes. This is used to calculate prop.bites */
      double max_sum = 0.0;
      for (int i = 0; i < Nnodes; i++) {
        double current_sum = u_cur[i * Nc + 3] + u_cur[i * Nc + 4] + u_cur[i * Nc + 5] + u_cur[i * Nc + 6] + u_cur[i * Nc + 10] + u_cur[i * Nc + 11] + u_cur[i * Nc + 12] + u_cur[i * Nc + 13] + u_cur[i * Nc + 17] + u_cur[i * Nc + 18] + u_cur[i * Nc + 19] + u_cur[i * Nc + 20] + u_cur[i * Nc + 22] + u_cur[i * Nc + 23] + u_cur[i * Nc + 24] + u_cur[i * Nc + 25];
        if (current_sum > max_sum) {
        max_sum = current_sum;
        }      
      } 
      
      v_new[12] = max_sum;
      v_new[13] = v_new[11]/v_new[12];


      /* Daily dispersal. 
       * Determine the pointer to the compartment state vector in the first node (u_cur). 
       * Use this to find the number of individuals at neighbours to the current node.
       * Clear the force of infection from neighboring nodes. Then iterate over all neighbors and update the value. 
       * Juvenile 'I_coupling' is the fourth item in the continuous state vector 'v'.
       * 'S_coupling' is the fifth item. Adult I_coupling and S_coupling are sixth and seventh. */
      v_new[4] = 0.0;
      v_new[5] = 0.0;
      v_new[6] = 0.0;
      v_new[7] = 0.0;
      v_new[14] = 0.0;
      v_new[15] = 0.0;
      v_new[16] = 0.0;
      v_new[17] = 0.0;

      /* Iterate over all neighbors. 
       * For breeding season, use breeding dispersal matrix,
       * during summer, use summer dispersal matrix (different ldata index). */
        
      if ((t<183) | ((t>456) & (t<548)) | ((t>821) & (t<913)) | ((t>1186) & (t<1278)) |
          ((t>1552) & (t<1644)) | ((t>1917) & (t<2009)) | ((t>2282) & (t<2374)) ) {

          for (int i = 0; i < Nnodes; i++) {

            double coupling = ldata[i+2588];

          if (u_cur[i * Nc + 5] > 0 & (u_cur[i * Nc + 3] + u_cur[i * Nc + 4] + u_cur[i * Nc + 5] + u_cur[i * Nc + 6]) > 0) {
            v_new[4] = v_new[4] + coupling * u_cur[i * Nc + 5]/(u_cur[i * Nc + 3] + u_cur[i * Nc + 4] + u_cur[i * Nc + 5] + u_cur[i * Nc + 6]);
          }

          if (u_cur[i * Nc + 2] > 0 & (u_cur[i * Nc + 3] + u_cur[i * Nc + 4] + u_cur[i * Nc + 5] + u_cur[i * Nc + 6]) > 0) {
            v_new[5] = v_new[5] + coupling * u_cur[i * Nc + 2]/(u_cur[i * Nc + 3] + u_cur[i * Nc + 4] + u_cur[i * Nc + 5] + u_cur[i * Nc + 6]);
          }

          if (u_cur[i * Nc + 12] > 0 & (u_cur[i * Nc + 10] + u_cur[i * Nc + 11] + u_cur[i * Nc + 12] + u_cur[i * Nc + 13]) > 0) {
            v_new[6] = v_new[6] + coupling * u_cur[i * Nc + 12]/(u_cur[i * Nc + 10] + u_cur[i * Nc + 11] + u_cur[i * Nc + 12] + u_cur[i * Nc + 13]);
          }

          if (u_cur[i * Nc + 2] > 0 & (u_cur[i * Nc + 10] + u_cur[i * Nc + 11] + u_cur[i * Nc + 12] + u_cur[i * Nc + 13]) > 0) {
            v_new[7] = v_new[7] + coupling * u_cur[i * Nc + 2]/(u_cur[i * Nc + 10] + u_cur[i * Nc + 11] + u_cur[i * Nc + 12] + u_cur[i * Nc + 13]);
          }

        }

    } else {

        for (int i = 0; i < Nnodes; i++) {

          double coupling = ldata[i+3986];

        if (u_cur[i * Nc + 5] > 0 & (u_cur[i * Nc + 3] + u_cur[i * Nc + 4] + u_cur[i * Nc + 5] + u_cur[i * Nc + 6]) > 0) {
          v_new[4] = v_new[4] + coupling * u_cur[i * Nc + 5]/(u_cur[i * Nc + 3] + u_cur[i * Nc + 4] + u_cur[i * Nc + 5] + u_cur[i * Nc + 6]);
        }

        if (u_cur[i * Nc + 2] > 0 & (u_cur[i * Nc + 3] + u_cur[i * Nc + 4] + u_cur[i * Nc + 5] + u_cur[i * Nc + 6]) > 0) {
          v_new[5] = v_new[5] + coupling * u_cur[i * Nc + 2]/(u_cur[i * Nc + 3] + u_cur[i * Nc + 4] + u_cur[i * Nc + 5] + u_cur[i * Nc + 6]);
        }

        if (u_cur[i * Nc + 12] > 0 & (u_cur[i * Nc + 10] + u_cur[i * Nc + 11] + u_cur[i * Nc + 12] + u_cur[i * Nc + 13]) > 0) {
          v_new[6] = v_new[6] + coupling * u_cur[i * Nc + 12]/(u_cur[i * Nc + 10] + u_cur[i * Nc + 11] + u_cur[i * Nc + 12] + u_cur[i * Nc + 13]);
        }

        if (u_cur[i * Nc + 2] > 0 & (u_cur[i * Nc + 10] + u_cur[i * Nc + 11] + u_cur[i * Nc + 12] + u_cur[i * Nc + 13]) > 0) {
          v_new[7] = v_new[7] + coupling * u_cur[i * Nc + 2]/(u_cur[i * Nc + 10] + u_cur[i * Nc + 11] + u_cur[i * Nc + 12] + u_cur[i * Nc + 13]);
        }
         
    }
  }

      /* Same principle for reservoir host, using different ldata. */

      if ((t<183) | ((t>456) & (t<548)) | ((t>821) & (t<913)) | ((t>1186) & (t<1278)) |
          ((t>1552) & (t<1644)) | ((t>1917) & (t<2009)) | ((t>2282) & (t<2374)) ) {

          for (int i = 0; i < Nnodes; i++) {

            double coupling = ldata[i+5384];

          if (u_cur[i * Nc + 19] > 0 & (u_cur[i * Nc + 17] + u_cur[i * Nc + 18] + u_cur[i * Nc + 19] + u_cur[i * Nc + 20]) > 0) {
            v_new[14] = v_new[14] + coupling * u_cur[i * Nc + 19]/(u_cur[i * Nc + 17] + u_cur[i * Nc + 18] + u_cur[i * Nc + 19] + u_cur[i * Nc + 20]);
          }

          if (u_cur[i * Nc + 2] > 0 & (u_cur[i * Nc + 17] + u_cur[i * Nc + 18] + u_cur[i * Nc + 19] + u_cur[i * Nc + 20]) > 0) {
            v_new[15] = v_new[15] + coupling * u_cur[i * Nc + 2]/(u_cur[i * Nc + 17] + u_cur[i * Nc + 18] + u_cur[i * Nc + 19] + u_cur[i * Nc + 20]);
          }

          if (u_cur[i * Nc + 24] > 0 & (u_cur[i * Nc + 22] + u_cur[i * Nc + 23] + u_cur[i * Nc + 24] + u_cur[i * Nc + 25]) > 0) {
            v_new[16] = v_new[16] + coupling * u_cur[i * Nc + 24]/(u_cur[i * Nc + 22] + u_cur[i * Nc + 23] + u_cur[i * Nc + 24] + u_cur[i * Nc + 25]);
          }

          if (u_cur[i * Nc + 2] > 0 & (u_cur[i * Nc + 22] + u_cur[i * Nc + 23] + u_cur[i * Nc + 24] + u_cur[i * Nc + 25]) > 0) {
            v_new[17] = v_new[17] + coupling * u_cur[i * Nc + 2]/(u_cur[i * Nc + 22] + u_cur[i * Nc + 23] + u_cur[i * Nc + 24] + u_cur[i * Nc + 25]);
          }


        }

    } else {

        for (int i = 0; i < Nnodes; i++) {

          double coupling = ldata[i+6782];

        if (u_cur[i * Nc + 19] > 0 & (u_cur[i * Nc + 17] + u_cur[i * Nc + 18] + u_cur[i * Nc + 19] + u_cur[i * Nc + 20]) > 0) {
            v_new[14] = v_new[14] + coupling * u_cur[i * Nc + 19]/(u_cur[i * Nc + 17] + u_cur[i * Nc + 18] + u_cur[i * Nc + 19] + u_cur[i * Nc + 20]);
        }

        if (u_cur[i * Nc + 2] > 0 & (u_cur[i * Nc + 17] + u_cur[i * Nc + 18] + u_cur[i * Nc + 19] + u_cur[i * Nc + 20]) > 0) {
            v_new[15] = v_new[15] + coupling * u_cur[i * Nc + 2]/(u_cur[i * Nc + 17] + u_cur[i * Nc + 18] + u_cur[i * Nc + 19] + u_cur[i * Nc + 20]);
        }

        if (u_cur[i * Nc + 24] > 0 & (u_cur[i * Nc + 22] + u_cur[i * Nc + 23] + u_cur[i * Nc + 24] + u_cur[i * Nc + 25]) > 0) {
            v_new[16] = v_new[16] + coupling * u_cur[i * Nc + 24]/(u_cur[i * Nc + 22] + u_cur[i * Nc + 23] + u_cur[i * Nc + 24] + u_cur[i * Nc + 25]);
        }
        if (u_cur[i * Nc + 2] > 0 & (u_cur[i * Nc + 22] + u_cur[i * Nc + 23] + u_cur[i * Nc + 24] + u_cur[i * Nc + 25]) > 0) {
            v_new[17] = v_new[17] + coupling * u_cur[i * Nc + 2]/(u_cur[i * Nc + 22] + u_cur[i * Nc + 23] + u_cur[i * Nc + 24] + u_cur[i * Nc + 25]);
        }

    }
  }
      
      /* Error check the new I_coupling value. Finally, if I_coupling
       * has changed compared to the previous value, return 1 to
       * indicate to the numerical solver that the transition rates must
       * be updated. */

      if (v_new[5] < 0.0)
          return SIMINF_ERR_V_IS_NEGATIVE;
      return 1; /* 1 if needs update */
  "

