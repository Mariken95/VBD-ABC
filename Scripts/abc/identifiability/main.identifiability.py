#!/usr/bin/python
import math
import sys
import os
import re
import time
#======================================================================
# functions
#======================================================================
def write_submission_script(filename, args, scenario, datasetn):
    submit_file = open(filename,'w')
    submit_str = "#!/bin/bash\n\n"
    submit_str += "#SBATCH --comment=identifiability \n"
    submit_str += "#SBATCH --partition=fat_rome \n"
    submit_str += "#SBATCH --time=00-26:00:00 \n"
    submit_str += "#SBATCH --job-name=abc_id_%d.%d \n" % (scenario, datasetn)
    submit_str += "#SBATCH --output=abc_id-%d.%d.txt \n" % (scenario, datasetn)
    submit_str += "#SBATCH --error=error_abc_id-%d.%d.txt \n" % (scenario, datasetn)
    submit_str += "#SBATCH --nodes=1 \n"
    submit_str += "#SBATCH --cpus-per-task=16 \n"
    submit_str += "#SBATCH --mail-type=FAIL \n"
    submit_str += "#SBATCH --mail-user=mariken.dewit@wur.nl \n \n"

    submit_str += "cd /gpfs/work4/0/prjs1002/usutu_sim/Scripts/ \n"

    submit_str += "module load 2023 \n"
    submit_str += "module load R/4.3.2-gfbf-2023a \n"

    submit_str += "srun R CMD BATCH --no-save --no-restore '--args%s' abc_run_slurm.R scen.%d.%d.out \n\n"  % (args, scenario, datasetn)

    submit_str += "exit 0 "

    submit_file.write(submit_str)
    submit_file.close()

#======================================================================
# Generate submission files
#======================================================================
parameterSet = [
    
# set A chain 1 
' userParameterFile=paste("abc_userParameters_id.A1.R"); previous_gens=NA; previous_epsilons=NA;',
# set B chain 1 
' userParameterFile=paste("abc_userParameters_id.B1.R"); previous_gens=NA; previous_epsilons=NA;', 
# set C chain 1 
' userParameterFile=paste("abc_userParameters_id.C1.R"); previous_gens=NA; previous_epsilons=NA;', 
# set D chain 1 
' userParameterFile=paste("abc_userParameters_id.D1.R"); previous_gens=NA; previous_epsilons=NA;', 
# set E chain 1 
' userParameterFile=paste("abc_userParameters_id.E1.R"); previous_gens=NA; previous_epsilons=NA;', 

# set A chain 2 
' userParameterFile=paste("abc_userParameters_id.A2.R"); previous_gens=NA; previous_epsilons=NA;',
# set B chain 2 
' userParameterFile=paste("abc_userParameters_id.B2.R"); previous_gens=NA; previous_epsilons=NA;', 
# set C chain 2 
' userParameterFile=paste("abc_userParameters_id.C2.R"); previous_gens=NA; previous_epsilons=NA;', 
# set D chain 2 
' userParameterFile=paste("abc_userParameters_id.D2.R"); previous_gens=NA; previous_epsilons=NA;', 
# set E chain 2 
' userParameterFile=paste("abc_userParameters_id.E2.R"); previous_gens=NA; previous_epsilons=NA;', 
]

for i in range(8,9):
    for j in range(0,1):
        submit_file_name = 'abc_ident_%d.%d.sh' % (i+1,j+1)
        scenario = parameterSet[i] + ' indexnr=%d;' % (j+1)
        write_submission_script(submit_file_name,scenario,i+1,j+1)
        print('submitting parameter set i:%d indexnr:%d' % (i+1,j+1))
        os.system("chmod +x %s" % (submit_file_name))
        os.system("sbatch %s" % (submit_file_name))
        time.sleep(.1)
