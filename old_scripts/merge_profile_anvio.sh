#!/bin/sh 
### Script to clean up the raw reads (trimming and host filtering)
### 2023-09-01 - Albert
### General options 
### -- set the job Name AND the job array --
#BSUB -J Merge
### -- specify queue -- 
#BSUB -q hpc
### -- ask for number of cores (default: 1) -- 
#BSUB -n 8
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 4:00  
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=4GB]"
### -- set the email address -- 
#BSUB -u apca@biosustain.dtu.dk
### -- send notification at completion -- 
#BSUB -N 
### -- Specify the output and error file. %J is the job-id %I is the job-array index -- 
#BSUB -o /work3/apca/orange_peel/04_anvio/log/Output_%J_%I.out 
#BSUB -e /work3/apca/orange_peel/04_anvio/log/Output_%J_%I.err 

### exit if any errors or unset variables are encountered
set +eu

### setting the working and output directory
PROJECT_DIR=/work3/apca/orange_peel
OUTPUT_DIR=$PROJECT_DIR/04_anvio
DB_DIR=$PROJECT_DIR/04_anvio/db

echo `date -Is` "-Start-"

# activate the Anvio environment
source /zhome/f5/1/58775/anaconda3/bin/activate anvio-7.1

### exit if any errors or unset variables are encountered
set -euo pipefail

# profile the samples to the enriched contigs database
# had to run the KEGG setup before running this: anvi-setup-kegg-kofams
#/zhome/f5/1/58775/anaconda3/envs/anvio-7.1/bin/anvi-merge ${OUTPUT_DIR}/profile/*/PROFILE.db \
#	-c ${DB_DIR}/contigs.db \
#	-o ${OUTPUT_DIR}/samples_merged

#/zhome/f5/1/58775/anaconda3/envs/anvio-7.1/bin/anvi-interactive -p ${OUTPUT_DIR}/samples_merged/PROFILE.db \
#	-c ${DB_DIR}/contigs.db

/zhome/f5/1/58775/anaconda3/envs/anvio-7.1/bin/anvi-estimate-metabolism -c ${DB_DIR}/contigs.db \
	--metagenome-mode

echo `date -Is` "-Finished-"