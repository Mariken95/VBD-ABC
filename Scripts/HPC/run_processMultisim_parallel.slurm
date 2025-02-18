#!/bin/bash
#SBATCH --comment=process_Multisim
#SBATCH --job-name="process_Multisim" 
#SBATCH --output=process_Multisim-%j.txt 
#SBATCH --error=error_process_Multisim-%j.txt 
#SBATCH --partition=fat_rome
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=10

cd /gpfs/work4/0/prjs1002/usutu_sim/Scripts

module load 2023 
module load R/4.3.2-gfbf-2023a
  
# Execute a series of serial programs (i.e. a 'pipeline'). This example pipeline consists of preprocessing, analysis and postprocessing.
# In this example, each of our programs take two arguments: the input (defined with -i argument) and an output (defined with -o).
# The '(' and ')' brackets are used to define a code block.
# The '&' sign is used to start each code block in the background, so that the pipeline starts running on input_file_1 to input_file_10 concurrently.
# Within the code block, the commands are executed sequentially. That way, the analysis will not start before the preprocessing is finished.
 
 
for i in `seq 1 $SLURM_NTASKS`; do
(
  Rscript output_processMultisim_parallel.R -i ../Output/Simulations/multi_20240805_A/multirun_$i.Rdata -o ../Output/Simulations/multi_20240805_A/processedRuns_$i
) &
done
wait