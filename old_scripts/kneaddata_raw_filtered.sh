#!/bin/sh 
### Script to clean up the raw reads (trimming and host filtering)
### 2023-08-10 - Albert
### General options 
### -- set the job Name AND the job array --
#BSUB -J My_array[1-4]
### -- specify queue -- 
#BSUB -q hpc
### -- ask for number of cores (default: 1) -- 
#BSUB -n 20
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
#BSUB -o /work3/apca/orange_peel/01_qc/log/Output_%J_%I.out 
#BSUB -e /work3/apca/orange_peel/01_qc/log/Output_%J_%I.err 

### exit if any errors or unset variables are encountered
set -euo pipefail

### setting the working and output directory
PROJECT_DIR=/work3/apca/orange_peel
OUTPUT_DIR=/work3/apca/orange_peel/01_qc

### Getting sample names and fq1 and fq2 files
id=${LSB_JOBINDEX}
SAMPLE=`cat sampletable | head -n $id | tail -1 | cut -f1`
FQ1=`cat sampletable | head -n $id | tail -1 | cut -f2`
FQ2=`cat sampletable | head -n $id | tail -1 | cut -f3`

### Creating needed output folders
if [ ! -d $OUTPUT_DIR/filtered ]; then
        mkdir $OUTPUT_DIR/filtered
fi

if [ ! -d $OUTPUT_DIR/filtered/$SAMPLE ]; then
        mkdir $OUTPUT_DIR/filtered/$SAMPLE
fi

echo `date -Is` "-Start-"

### Filtering the raw reads by trimming low quality reads and filtering out host reads
cd /zhome/f5/1/58775/.local/bin/
kneaddata -i1 ${FQ1} -i2 ${FQ2} \
	--reference-db ${PROJECT_DIR}/orange_genome \
	--output ${OUTPUT_DIR}/filtered/${SAMPLE} \
	--bypass-trim \
	--trf /zhome/f5/1/58775/anaconda3/bin
	
#	--bypass-trim \
#	--bypass-trf \
#	--bowtie2 /zhome/f5/1/58775/.local/bin
#	--trimmomatic /zhome/f5/1/58775/.local/bin/ \
#	--trimmomatic-options="MINLEN:60 ILLUMINACLIP:/zhome/f5/1/58775/.local/lib/python2.7/site-packages/kneaddata/adapters/TruSeq3-PE.fa:2:30:10:8:TRUE SLIDINGWINDOW:4:20 MINLEN:75" \
#	--bowtie2-path /zhome/f5/1/58775/.local/bin
#	--trf /zhome/f5/1/58775/.local/bin/
	

echo `date -Is` "-Finished-"