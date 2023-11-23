#!/bin/sh 
### Script to determine how much Prokaryote and Eukaryote we have in the samples 
### 2023-08-28 - Albert
### General options 
### -- set the job Name AND the job array --
#BSUB -J WhoK_array[1-4]
### -- specify queue -- 
#BSUB -q hpc
### -- ask for number of cores (default: 1) -- 
#BSUB -n 4
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
#BSUB -o /work3/apca/orange_peel/05_whokaryote/log/Output_%J_%I.out 
#BSUB -e /work3/apca/orange_peel/05_whokaryote/log/Output_%J_%I.err 

### exit if any errors or unset variables are encountered
set +eu

### setting the working and output directory
PROJECT_DIR=/work3/apca/orange_peel
OUTPUT_DIR=/work3/apca/orange_peel/05_whokaryote

### Getting sample names and fq1 and fq2 files
id=${LSB_JOBINDEX}
SAMPLE=`cat sampletable | head -n $id | tail -1 | cut -f1`
CONTIG=$PROJECT_DIR/02_assembly/assembly/$SAMPLE/final.contigs.fa


### Creating needed output folders
if [ ! -d $OUTPUT_DIR/whokaryote ]; then
        mkdir $OUTPUT_DIR/whokaryote
fi

if [ ! -d $OUTPUT_DIR/whokaryote/$SAMPLE ]; then
        mkdir $OUTPUT_DIR/whokaryote/$SAMPLE
fi

echo `date -Is` "-Start-"

# need to activate whokaryote environment
source /zhome/f5/1/58775/anaconda3/bin/activate whokaryote

set -euo pipefail

### Running the whokaryote
/zhome/f5/1/58775/anaconda3/envs/whokaryote/bin/whokaryote.py --contigs ${CONTIG} \
	--minsize 1000 \
	--outdir ${OUTPUT_DIR}/whokaryote/${SAMPLE}

echo `date -Is` "-Finished-"