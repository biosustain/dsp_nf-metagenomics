#!/bin/sh 
### Script to analyse the quality of the sequenced data before doing further quality controls or analysis
### 2023-08-15 - Albert
### General options 
### -- set the job Name AND the job array --
#BSUB -J My_array[1-4]
### -- specify queue -- 
#BSUB -q hpc
### -- ask for number of cores (default: 1) -- 
#BSUB -n 1
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
#BSUB -o /work3/apca/orange_peel/00_fastqc/log/Output_%J_%I.out 
#BSUB -e /work3/apca/orange_peel/00_fastqc/log/Output_%J_%I.err 

### exit if any errors or unset variables are encountered
set -euo pipefail

### setting the working and output directory
PROJECT_DIR=/work3/apca/orange_peel
OUTPUT_DIR=/work3/apca/orange_peel/00_fastqc

### Getting sample names and fq1 and fq2 files
id=${LSB_JOBINDEX}
SAMPLE=`cat sampletable | head -n $id | tail -1 | cut -f1`
FQ1=`cat sampletable | head -n $id | tail -1 | cut -f2`
FQ2=`cat sampletable | head -n $id | tail -1 | cut -f3`

### Creating needed output folders
if [ ! -d $OUTPUT_DIR/fastqc ]; then
        mkdir $OUTPUT_DIR/fastqc
fi

if [ ! -d $OUTPUT_DIR/fastqc/$SAMPLE ]; then
        mkdir $OUTPUT_DIR/fastqc/$SAMPLE
fi

echo `date -Is` "-Start-"

### Creating fastQC reports for each fastq file
/zhome/f5/1/58775/anaconda3/bin/fastqc -o $OUTPUT_DIR/fastqc/$SAMPLE ${FQ1} ${FQ2}
# Would be good to remove zip files after or mv the files somewhere whenever they were checked
# rm  $OUTPUT_DIR/fastqc/$SAMPLE/*.zip

echo `date -Is` "-Finished-"