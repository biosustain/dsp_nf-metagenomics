#!/bin/sh 
### Script to determine how much Prokaryote and Eukaryote we have in the samples 
### 2023-08-29 - Albert
### General options 
### -- set the job Name AND the job array --
#BSUB -J ReformatContigs
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
#BSUB -o /work3/apca/orange_peel/03_reformat_check_contigs/log/Output_%J_%I.out 
#BSUB -e /work3/apca/orange_peel/03_reformat_check_contigs/log/Output_%J_%I.err 

### exit if any errors or unset variables are encountered
set +eu

### setting the working and output directory
PROJECT_DIR=/work3/apca/orange_peel
OUTPUT_DIR=/work3/apca/orange_peel/03_reformat_check_contigs

### Building the contig variable
CONTIG=$PROJECT_DIR/02_assembly/coassembly/final.contigs.fa

### Creating needed reformatted contigs folders
if [ ! -d $OUTPUT_DIR/contigs ]; then
        mkdir $OUTPUT_DIR/contigs
fi

echo `date -Is` "-Start-"

# activate the Anvio environment
source /zhome/f5/1/58775/anaconda3/bin/activate anvio-7.1

### exit if any errors or unset variables are encountered
set -euo pipefail


# contigs shorter than 2500 bp are removed as they are too short according to anvio tutorial
/zhome/f5/1/58775/anaconda3/envs/anvio-7.1/bin/anvi-script-reformat-fasta ${CONTIG} \
	-o ${OUTPUT_DIR}/contigs/contigs.fa \
	-l 2500 \
	--simplify-names \
	-r ${OUTPUT_DIR}/name_conversions.txt

echo `date -Is` "-Finished-"