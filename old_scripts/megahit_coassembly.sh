#!/bin/sh 
### Script to coassemble all samples in preparation to run Anvio
### 2023-08-17 - Albert
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J Coassembly
### -- ask for number of cores (default: 1) -- 
#BSUB -n 28 
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
#BSUB -o /work3/apca/orange_peel/02_assembly/log/Output_%J.out 
#BSUB -e /work3/apca/orange_peel/02_assembly/log/Output_%J.err 


### exit if any errors or unset variables are encountered
set -euo pipefail

### setting the working and output directory
PROJECT_DIR=/work3/apca/orange_peel
OUTPUT_DIR=/work3/apca/orange_peel/02_assembly

### Getting sample names and fq1 and fq2 files
R1s=$PROJECT_DIR/01_qc/filtered/O1/O1_1_kneaddata_paired_1.fastq,$PROJECT_DIR/01_qc/filtered/O2/O2_1_kneaddata_paired_1.fastq,$PROJECT_DIR/01_qc/filtered/O3/O3_1_kneaddata_paired_1.fastq,$PROJECT_DIR/01_qc/filtered/O4/O4_1_kneaddata_paired_1.fastq
R2s=$PROJECT_DIR/01_qc/filtered/O1/O1_1_kneaddata_paired_2.fastq,$PROJECT_DIR/01_qc/filtered/O2/O2_1_kneaddata_paired_2.fastq,$PROJECT_DIR/01_qc/filtered/O3/O3_1_kneaddata_paired_2.fastq,$PROJECT_DIR/01_qc/filtered/O4/O4_1_kneaddata_paired_2.fastq

### Creating needed output folders
if [ ! -d $OUTPUT_DIR/coassembly ]; then
        mkdir $OUTPUT_DIR/coassembly
fi

echo `date -Is` "-Start-"

### Building the orange genome database to map reads onto it
/zhome/f5/1/58775/anaconda3/bin/megahit -1 ${R1s} -2 ${R2s} \
	–min-contig-len 1000 \
	-m 0.85 \
	-o $OUTPUT_DIR/coassembly \
	-t 28


echo `date -Is` "-Finished-"