#!/bin/sh 
### Script to clean up the raw reads (trimming and host filtering)
### 2023-09-01 - Albert
### General options 
### -- set the job Name AND the job array --
#BSUB -J Profile[1-4]
### -- specify queue -- 
#BSUB -q hpc
### -- ask for number of cores (default: 1) -- 
#BSUB -n 32
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 24:00  
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

### Getting sample names and fq1 and fq2 files
id=${LSB_JOBINDEX}
SAMPLE=`cat sampletable | head -n $id | tail -1 | cut -f1`
BAM=${PROJECT_DIR}/03_mapping/mapping/${SAMPLE}/${SAMPLE}.bam

### Creating needed output folders
if [ ! -d $OUTPUT_DIR/profile ]; then
        mkdir $OUTPUT_DIR/profile
fi

echo `date -Is` "-Start-"

# activate the Anvio environment
source /zhome/f5/1/58775/anaconda3/bin/activate anvio-7.1

### exit if any errors or unset variables are encountered
set -euo pipefail

# profile the samples to the enriched contigs database
/zhome/f5/1/58775/anaconda3/envs/anvio-7.1/bin/anvi-profile -i ${PROJECT_DIR}/03_mapping/mapping/${SAMPLE}/${SAMPLE}.bam \
	-c ${DB_DIR}/contigs.db \
	-o ${OUTPUT_DIR}/profile/${SAMPLE}

echo `date -Is` "-Finished-"