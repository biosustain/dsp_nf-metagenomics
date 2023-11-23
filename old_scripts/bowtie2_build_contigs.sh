#!/bin/sh 
### Script to generate an orange genome database to map reads onto it
### 2023-08-09 - Albert
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J bowtie2_contigs
### -- ask for number of cores (default: 1) -- 
#BSUB -n 8 
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=4GB]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot -- 
#BSUB -M 4GB
### -- set walltime limit: hh:mm -- 
#BSUB -W 24:00 
### -- set the email address -- 
#BSUB -u apca@biosustain.dtu.dk
### -- send notification at start -- 
#BSUB -B 
### -- send notification at completion -- 
#BSUB -N 
### -- Specify the output and error file. %J is the job-id -- 
### -- -o and -e mean append, -oo and -eo mean overwrite -- 
#BSUB -o /work3/apca/orange_peel/03_mapping/log/Output_%J.out 
#BSUB -e /work3/apca/orange_peel/03_mapping/log/Output_%J.err 

### exit if any errors or unset variables are encountered
set -euo pipefail

### setting the working directory
PROJECT_DIR=/work3/apca/orange_peel
DB_DIR=/work3/apca/orange_peel/03_mapping

### Building the orange genome database to map reads onto it. This version of bowtie2 downloaded with kneaddata does not support multithread
/zhome/f5/1/58775/.local/bin/bowtie2-build $PROJECT_DIR/03_reformat_check_contigs/contigs/contigs.fa $DB_DIR/contigs 

echo `date -Is` "-Finished-"