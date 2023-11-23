#!/bin/sh 
### Script to clean up the raw reads (trimming and host filtering)
### 2023-08-11 - Albert
### General options 
### -- set the job Name AND the job array --
#BSUB -J MegaHit_array[1-4]
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
#BSUB -o /work3/apca/orange_peel/02_assembly/log/Output_%J_%I.out 
#BSUB -e /work3/apca/orange_peel/02_assembly/log/Output_%J_%I.err 

### exit if any errors or unset variables are encountered
set -euo pipefail

### setting the working and output directory
PROJECT_DIR=/work3/apca/orange_peel
OUTPUT_DIR=/work3/apca/orange_peel/02_assembly

### Getting sample names and fq1 and fq2 files
id=${LSB_JOBINDEX}
SAMPLE=`cat sampletable | head -n $id | tail -1 | cut -f1`
FQ1=$PROJECT_DIR/01_qc/filtered/$SAMPLE/${SAMPLE}_1_kneaddata_paired_1.fastq
FQ2=$PROJECT_DIR/01_qc/filtered/$SAMPLE/${SAMPLE}_1_kneaddata_paired_2.fastq

### Creating needed output folders
if [ ! -d $OUTPUT_DIR/assembly ]; then
        mkdir $OUTPUT_DIR/assembly
fi

#if [ ! -d $OUTPUT_DIR/assembly/$SAMPLE ]; then
#        mkdir $OUTPUT_DIR/assembly/$SAMPLE
#fi

echo `date -Is` "-Start-"

### Building the orange genome database to map reads onto it
/zhome/f5/1/58775/anaconda3/bin/megahit -1 ${FQ1} -2 ${FQ2} \
	-m 0.9 -t 8 -o ${OUTPUT_DIR}/assembly/${SAMPLE}

echo `date -Is` "-Finished-"