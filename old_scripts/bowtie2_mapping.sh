#!/bin/sh 
### Script to clean up the raw reads (trimming and host filtering)
### 2023-08-31 - Albert
### General options 
### -- set the job Name AND the job array --
#BSUB -J Map_array[1-4]
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
#BSUB -o /work3/apca/orange_peel/03_mapping/log/Output_%J_%I.out 
#BSUB -e /work3/apca/orange_peel/03_mapping/log/Output_%J_%I.err 

### exit if any errors or unset variables are encountered
set +eu

### setting the working and output directory
PROJECT_DIR=/work3/apca/orange_peel
OUTPUT_DIR=$PROJECT_DIR/03_mapping
DB_DIR=$PROJECT_DIR/03_mapping/contigs

### Getting sample names and fq1 and fq2 files
id=${LSB_JOBINDEX}
SAMPLE=`cat sampletable | head -n $id | tail -1 | cut -f1`
FQ1=$PROJECT_DIR/01_qc/filtered/$SAMPLE/${SAMPLE}_1_kneaddata_paired_1.fastq
FQ2=$PROJECT_DIR/01_qc/filtered/$SAMPLE/${SAMPLE}_1_kneaddata_paired_2.fastq

### Creating needed output folders
if [ ! -d $OUTPUT_DIR/mapping ]; then
        mkdir $OUTPUT_DIR/mapping
fi

if [ ! -d $OUTPUT_DIR/mapping/$SAMPLE ]; then
        mkdir $OUTPUT_DIR/mapping/$SAMPLE
fi

echo `date -Is` "-Start-"

# activate the Anvio environment
source /zhome/f5/1/58775/anaconda3/bin/activate anvio-7.1

### exit if any errors or unset variables are encountered
set -euo pipefail

# getting an indexed bam file for each of our samples
# first, creating a sam file
# second, compressing them to bam files
# third, transforming the raw bam files to a bam file that Anvio can understand
/zhome/f5/1/58775/anaconda3/envs/anvio-7.1/bin/bowtie2 --threads 8 \
	-x ${DB_DIR} \
	-1 ${FQ1} -2 ${FQ2} \
	-S ${OUTPUT_DIR}/mapping/${SAMPLE}/${SAMPLE}.sam &&
 
/zhome/f5/1/58775/anaconda3/envs/anvio-7.1/bin/samtools view -F 4 \
	-bS ${OUTPUT_DIR}/mapping/${SAMPLE}/${SAMPLE}.sam > ${OUTPUT_DIR}/mapping/${SAMPLE}/${SAMPLE}-RAW.bam &&

/zhome/f5/1/58775/anaconda3/envs/anvio-7.1/bin/anvi-init-bam ${OUTPUT_DIR}/mapping/${SAMPLE}/${SAMPLE}-RAW.bam \
	-o ${OUTPUT_DIR}/mapping/${SAMPLE}/${SAMPLE}.bam

echo `date -Is` "-Finished-"