#!/bin/sh 
### Script to clean up the raw reads (trimming and host filtering)
### 2023-08-31 - Albert
### General options 
### -- set the job Name AND the job array --
#BSUB -J EnrichDB
### -- specify queue -- 
#BSUB -q hpc
### -- ask for number of cores (default: 1) -- 
#BSUB -n 8
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

### Creating needed output folders
if [ ! -d $OUTPUT_DIR/db ]; then
        mkdir $OUTPUT_DIR/db
fi

echo `date -Is` "-Start-"

# activate the Anvio environment
source /zhome/f5/1/58775/anaconda3/bin/activate anvio-7.1

### exit if any errors or unset variables are encountered
set -euo pipefail

# contigs fasta file converted into an Anvio SQL database
# Then we enrich the database:
# first, annotation of HMMS
# second, annotation of COGs
# third, annotation of taxonomy
/zhome/f5/1/58775/anaconda3/envs/anvio-7.1/bin/anvi-gen-contigs-database -f ${PROJECT_DIR}/03_reformat_check_contigs/contigs/contigs.fa \
	-o ${DB_DIR}/contigs.db \
	-n "Orange Peel contigs database"

/zhome/f5/1/58775/anaconda3/envs/anvio-7.1/bin/anvi-run-hmms --just-do-it \
	-c ${DB_DIR}/contigs.db \
	--num-threads 8
 
/zhome/f5/1/58775/anaconda3/envs/anvio-7.1/bin/anvi-run-ncbi-cogs -c ${DB_DIR}/contigs.db \
	--num-threads 8

/zhome/f5/1/58775/anaconda3/envs/anvio-7.1/bin/anvi-run-scg-taxonomy -c ${DB_DIR}/contigs.db \
	--num-threads 8

/zhome/f5/1/58775/anaconda3/envs/anvio-7.1/bin/anvi-run-kegg-kofams --just-do-it \
	-c ${DB_DIR}/contigs.db \
	--num-threads 8

echo `date -Is` "-Finished-"