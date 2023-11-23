#!/bin/sh 
### Script to generate an orange genome database to map reads onto it
### 2023-08-30 - Albert
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J run_whok_stats
### -- ask for number of cores (default: 1) -- 
#BSUB -n 1 
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=4GB]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot -- 
#BSUB -M 4GB
### -- set walltime limit: hh:mm -- 
#BSUB -W 4:00 
### -- set the email address -- 
#BSUB -u apca@biosustain.dtu.dk
### -- send notification at start -- 
#BSUB -B 
### -- send notification at completion -- 
#BSUB -N 
### -- Specify the output and error file. %J is the job-id -- 
### -- -o and -e mean append, -oo and -eo mean overwrite -- 
#BSUB -o /work3/apca/orange_peel/05_whokaryote/log/Output_%J.out 
#BSUB -e /work3/apca/orange_peel/05_whokaryote/log/Output_%J.err 

### exit if any errors or unset variables are encountered
set -euo pipefail

### setting the working directory
PROJECT_DIR=/work3/apca/orange_peel

python ${PROJECT_DIR}/src/whokaryote_stats.py

# THIS DOES NOT WORK INSIDE THE SHELL SCRIPT BUT IT WORKS OUTSIDE SO I GO AHEAD FOR NOW