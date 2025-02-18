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
    submit_str += "#SBATCH --comment=multirun \n"
    submit_str += "#SBATCH --partition=fat_rome \n"
    submit_str += "#SBATCH --time=00-00:20:00 \n"
    submit_str += "#SBATCH --job-name=sim_%d.%d \n" % (scenario, datasetn)
    submit_str += "#SBATCH --output=sim-%d.%d.txt \n" % (scenario, datasetn)
    submit_str += "#SBATCH --error=error_sim-%d.%d.txt \n" % (scenario, datasetn)
    submit_str += "#SBATCH --mail-type=FAIL \n"
    submit_str += "#SBATCH --mail-user=mariken.dewit@wur.nl \n \n"

    submit_str += "cd /gpfs/work4/0/prjs1002/usutu_sim/Scripts/ \n"

    submit_str += "module load 2023 \n"
    submit_str += "module load R/4.3.2-gfbf-2023a \n"

    submit_str += "srun R CMD BATCH --no-save --no-restore '--args%s' main_run_multisim.R run.%d.%d.out \n\n"  % (args, scenario, datasetn)

    submit_str += "exit 0 "

    submit_file.write(submit_str)
    submit_file.close()

#======================================================================
# Generate submission files
#======================================================================
version = [
    # model A
' nr.runs=10; scenario.dispersal=paste("blackbird"); scenario.parameters=paste("A"); posterior=paste("posteriorA"); outputPath=paste("../Output/Simulations/multi_20240805_A/");',

    # model B
' nr.runs=10; scenario.dispersal=paste("sc20"); scenario.parameters=paste("B"); posterior=paste("posteriorB"); outputPath=paste("../Output/Simulations/multi_20240805_B/");',

    # model C
' nr.runs=10; scenario.dispersal=paste("blackbird"); scenario.parameters=paste("C"); posterior=paste("posteriorC"); outputPath=paste("../Output/Simulations/multi_20240805_C/");',

    # model D
' nr.runs=10; scenario.dispersal=paste("sc20"); scenario.parameters=paste("D"); posterior=paste("posteriorD"); outputPath=paste("../Output/Simulations/multi_20240805_D/");',

    # model S1
' nr.runs=10; scenario.dispersal=paste("sc20"); scenario.parameters=paste("S1"); posterior=paste("posteriorS1"); outputPath=paste("../Output/Simulations/multi_20240805_S1/");',

    # model S2
' nr.runs=10; scenario.dispersal=paste("sc20"); scenario.parameters=paste("S2"); posterior=paste("posteriorS2"); outputPath=paste("../Output/Simulations/multi_20240805_S2/");',

    # model S3
' nr.runs=10; scenario.dispersal=paste("sensitivity"); scenario.parameters=paste("S3"); posterior=paste("posteriorS3"); outputPath=paste("../Output/Simulations/multi_20240805_S3/");',

]

for i in range(2,3):
    for j in range(0,10):
        submit_file_name = 'multirun_%d.%d.sh' % (i+1,j+1)
        scenario = version[i] + ' indexnr=%d;' % (j+1)
        write_submission_script(submit_file_name,scenario,i+1,j+1)
        print('submitting multiple simulations i:%d indexnr:%d' % (i+1,j+1))
        os.system("chmod +x %s" % (submit_file_name))
        os.system("sbatch %s" % (submit_file_name))
        time.sleep(.1)
