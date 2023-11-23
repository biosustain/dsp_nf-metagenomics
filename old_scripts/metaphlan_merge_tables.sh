#!/bin/sh 
### Script to clean up the raw reads (trimming and host filtering)
### 2023-09-13 - Albert
### General options 
### -- set the job Name AND the job array --
#BSUB -J MergeTables
### -- specify queue -- 
#BSUB -q hpc
### -- ask for number of cores (default: 1) -- 
#BSUB -n 1
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 12:00  
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=30GB]"
### -- set the email address -- 
#BSUB -u apca@biosustain.dtu.dk
### -- send notification at completion -- 
#BSUB -N 
### -- Specify the output and error file. %J is the job-id %I is the job-array index -- 
#BSUB -o /work3/apca/metaphlan/02_merge_tables/log/Output_%J_%I.out 
#BSUB -e /work3/apca/metaphlan/02_merge_tables/log/Output_%J_%I.err 

### exit if any errors or unset variables are encountered
set +eu

PROFILES="/work3/apca/metaphlan/01_metaphlan/output/"
OUT="/work3/apca/metaphlan/02_merge_tables/"

echo `date` "Starting!"

# activate the Anvio environment
source /zhome/f5/1/58775/anaconda3/bin/activate mpa

### exit if any errors or unset variables are encountered
set -euo pipefail

/zhome/f5/1/58775/anaconda3/envs/mpa/bin/merge_metaphlan_tables.py ${PROFILES}/*profiled_metagenome.txt > ${OUT}/merged_abundance_table.txt

echo `date` Finished!
