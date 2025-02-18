# VBD-ABC
This repository contains the code associated with the paper:

Mariken M. de Wit, Gaël Beaunée*, Martha Dellar*, Emmanuelle Münger*, Louie Krol, Nnomzie Atama, Jurrian van Irsel, Henk van der Jeugd, Judith M.A. van den Brand, Chantal Reusken, Marion Koopmans, Mart C.M. de Jong, Reina S. Sikkema, Quirine ten Bosch - **"Silent reservoirs are shaping disease emergence: 
the case of Usutu virus in the Netherlands"**. The paper has been published in xx.

All code contained within this repository is released under the xxx licence.

# Overview
This code creates spatio-temporal transmission simulation models for vector-borne diseases. These models are fitted to observational data using Approximate Bayesian Computation. After fitting of the model to the data, simulations are performed using posterior estimates and several outbreak-related metrics are calculated to reconstruct the outbreak. The simulation models are built using the SimInf model framework and applied to Usutu virus in the Netherlands.

# System requirements
Data preparation and output analyses require a standard computer. Running of simulations and fitting to data requires a high-performance cluster (HPC) system.
The software has been tested on 1) standard computer: Windows 10, R 4.2.1 and 2) HPC Snellius system using fcn nodes provided by SURF: Linux, R 4.3.2. 

Users should install the following packages before executing the code:

```
install.packages(c('tidyverse', 'lubridate', 'patchwork', 'scales','SimInf', 'reshape2', 'ggridges', 'RColorBrewer', 'HDInterval', 'sf', 'raster', 'stars', 'readxl', 'tmap', 'spatialEco', 'sp', 'tmap', 'data.table', 'matlib')
```
The packages should install in a few minutes on a standard computer. 

The versions of packages the code has been tested on:
```
matlib_0.9.6       data.table_1.14.10 spatialEco_2.0-2   tmap_3.3-4         readxl_1.4.3       stars_0.6-4        abind_1.4-5        raster_3.6-26      sp_2.1-2           sf_1.0-14          HDInterval_0.2.4   RColorBrewer_1.1-3 ggridges_0.5.4     reshape2_1.4.4     SimInf_9.8.1       scales_1.3.0       patchwork_1.3.0    lubridate_1.9.3    forcats_1.0.0      stringr_1.5.1      dplyr_1.1.2        purrr_1.0.2        readr_2.1.4        tidyr_1.3.0        tibble_3.2.1       ggplot2_3.5.1      tidyverse_2.0.0
```  

Running simulations and fitting the models on the HPC requires the installation of the following packages:
```
install.packages(c('tidyverse', SimInf', 'zoo', 'data.table', 'optparse', 'lhs', 'scales', 'progress', 'processx', 'filelock','philentropy','lubridate'))
``` 

# Instructions for use

The folder `Scripts` contains code to run analyses and reproduce figures in the manuscript and supplements. Files starting with 'abc_' relate to the model calibration steps. Files starting with 'prep_' contain data preparation steps for model inputs and have associated 'function_' files with the same extension when relevant. 'main_' files are used to run the models (using associated function files) & create plots.

The folder `Scripts\HPC` contains run scripts to execute runs on the HPC. The folder `Scripts\abc` contains scripts related to the ABC algorithm for model fitting.

The folder `Data` contains input data used in Scripts.

## Data preparation
The model raster is created from _prep_createRaster.R_. Bird abundance, mosquito abundance, and temperature data are prepared in _prep_inputData_ with the associated function file _functions_inputData.R_.
These input data are further transformed into model input in _prep_processData.R_ (& _functions_processData.R_) resulting in starting conditions, events dataframes, and bird movement probabilities.
The files _prep_parameters_x.R_ contain all input parameters for the models (there are versions for both single-host and two-host models). This files contains parameters, compartments, equations, and the SimInf-related objects: E matrix, N matrix, pts_fun, v0.

Visualisations of bird dispersal input as shown in the supplementary material can be created from _supp_birdDispersalSimulation.R_.

### Input data
Temperature data can be obtained from the KNMI data platform: https://dataplatform.knmi.nl/dataset/tg1-5.

Bird abundance data can be obtained from _Dellar M, Sierdsema H, van Bodegom PM, Schrama MJJ, Geerling G (2024) The future abundance of key bird species for pathogen transmission in the Netherlands. bioRxiv 2024.10.14.617534_.

Mosquito abundance data can be obtained from _Krol L, Dellar M, Ibáñez-Justicia A, Schrier G van der, Schrama M, Geerling GW, Bodegom PM van (2024) Combined effects of future climate and land use change on mosquitoes: The distribution of Culex pipiens under One Health scenarios in the Netherlands. https://doi.org/10.21203/RS.3.RS-5298493/V1_

Mosquito abundance data are available upon request to the authors. To ensure reproducibility of our results, we therefore added the processed version of these data so the model can still be run without the original dataset. These can be found in `data/model_input/startingabundance.RData` & `data/model_input/events_5kmNoWater_2016_22_agestruc_dispersal.bb.RData`)


## Running the model
The model is run for a single iteration through the _main_run_model_ script using _functions_run_model_x.R_. This is set up to run on a local machine and will take a few minutes to run. Multiple iterations can be run through _main_run_multisim.R_ which is set up to be run on a Linux high-performance cluster (HPC). This script uses the same function file as the 'single iteration' script. These runs are submitted to the HPC through _launch.multirun.py_.

## Analysing output
Results from multiple model iterations can be processed by running _output_processMultirun_parallel.R_ on an HPC (by using _run_processMultirun_parallel.slurm_), after which plots can be generated from _main_plots.R_. Functions for the processing steps can be found in _functions_processOutput.R_.

## Model calibration
All files related to model calibration start with 'abc_'.  A more detailed explanation regarding the implementation of the algorithm can be found on https://gaelbn.github.io/BRREWABC/articles/BRREWABC.html. ABC calibration runs can be started from _main.abc.py_. The algorithm requires 5 files that are the same for all model versions: 'run_slurm.R', 'subjob.R', 'function_parallel.R', 'defaultParameters.R', 'updateParameters.R. Additionally, there is a 'defineModel_x.R' script for each model version and an 'userParameters_x.R' script for each model version and each chain. 

Summary statistics are calculated from simulation output using _abc_functions_summaryStatistics.R_. Calibration output is analysed in _abc_diagnostics.R_

The model calibration was first validated by parameter values from simulated data. Simulated datasets were created in _abc.validation_create.ssobs.R_.
The same abc-files (as described above) are used for fitting  to simulated data. These scripts can be found in `Scripts/abc/identifiability` and the runs can be launched from _main.identifiability.py_ in this folder.