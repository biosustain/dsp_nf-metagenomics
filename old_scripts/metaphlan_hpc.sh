#!/bin/sh 
### Script to clean up the raw reads (trimming and host filtering)
### 2023-08-17 - Albert
### General options 
### -- set the job Name AND the job array --
#BSUB -J Metaphlan_array[1-4]
### -- specify queue -- 
#BSUB -q hpc
### -- ask for number of cores (default: 1) -- 
#BSUB -n 8
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- set walltime limit: hh:mm -- 
#BSUB -W 12:00  
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=4GB]"
### -- set the email address -- 
#BSUB -u apca@biosustain.dtu.dk
### -- send notification at completion -- 
#BSUB -N 
### -- Specify the output and error file. %J is the job-id %I is the job-array index -- 
#BSUB -o /work3/apca/metaphlan/01_metaphlan/log/Output_%J_%I.out 
#BSUB -e /work3/apca/metaphlan/01_metaphlan/log/Output_%J_%I.err 

### exit if any errors or unset variables are encountered
set +eu

# to get a sample table with the path to the high quality filtered reads
# sed 's/samples/01_qc\/filtered/g' sampletable > sampletable_hiqc
# sed -i 's/1\.fq\.gz/1_kneaddata_paired_1\.fastq/g' sampletable_hiqc
# sed -i 's/2\.fq\.gz/1_kneaddata_paired_2\.fastq/g' sampletable_hiqc
readonly SAMPLETABLE="sampletable_hiqc"
readonly UNIQUESAMPLENAMES="samplelist"
readonly BOWTIEDIR="01_metaphlan/bowtie"
readonly OUTDIR="01_metaphlan/output"

mkdir -p "${BOWTIEDIR}"  "${OUTDIR}"

## read from a single line corresponding to the job array ID to get the sample name to process
#read SAMPLEID < (cat ${UNIQUESAMPLENAMES} | sed "${LSB_JOBINDEX}q;d")
SAMPLEID=`cat ${SAMPLETABLE} | head -n ${LSB_JOBINDEX} | tail -1 | cut -f1`

## identify all relevant sampletable lines, extract fastq paths
FQ1=$(grep "^${SAMPLEID}[[:space:]]" ${SAMPLETABLE} | tr " " "\t" | cut -f 2 | paste -s -d ",")
FQ2=$(grep "^${SAMPLEID}[[:space:]]" ${SAMPLETABLE} | tr " " "\t" | cut -f 3 | paste -s -d ",")

BOWTIEFILE="${BOWTIEDIR}/${SAMPLEID}.bowtie2.bz2"
OUTFILE="${OUTDIR}/${SAMPLEID}.profiled_metagenome.txt"
TEMPFILE="${OUTDIR}/${SAMPLEID}.incomplete.txt"
OUTFILE2="${OUTDIR}/${SAMPLEID}.marker_counts.txt"
OUTFILE3="${OUTDIR}/${SAMPLEID}.reads_map.txt"
OUTFILE4="${OUTDIR}/${SAMPLEID}.marker_ab_table.txt"

echo `date` "Starting ${SAMPLEID}"

# activate the Anvio environment
source /zhome/f5/1/58775/anaconda3/bin/activate mpa

### exit if any errors or unset variables are encountered
set -euo pipefail

echo "$FQ1"
echo "$FQ2"

[ -f "${OUTFILE}" ] && { echo "Stopping because output file already exists" ; exit 0; } 

## Remove previously created temp and bowtie files, if they exist
rm -f  "${TEMPFILE}"  "${BOWTIEFILE}"

/zhome/f5/1/58775/anaconda3/envs/mpa/bin/metaphlan "${FQ1},${FQ2}" \
	-t rel_ab \
	--bowtie2db /work3/apca/metaphlan/databases \
    --bowtie2out "${BOWTIEFILE}" \
    --nproc 8 \
    --input_type fastq \
    > "${TEMPFILE}"

mv "${TEMPFILE}"  "${OUTFILE}"

echo `date` adding marker_counts, reads_map, and marker_ab_table

/zhome/f5/1/58775/anaconda3/envs/mpa/bin/metaphlan "${BOWTIEFILE}" \
	-t marker_counts \
	--bowtie2db /work3/apca/metaphlan/databases \
    --nproc 8 \
    --input_type bowtie2out \
    > "${TEMPFILE}"


## only keep marker genes with hits
sed -i '/\t0$/d' "${TEMPFILE}"

mv "${TEMPFILE}"  "${OUTFILE2}"

/zhome/f5/1/58775/anaconda3/envs/mpa/bin/metaphlan "${BOWTIEFILE}" \
	-t reads_map \
	--bowtie2db /work3/apca/metaphlan/databases \
    --nproc 8 \
    --input_type bowtie2out \
    > "${TEMPFILE}"

mv "${TEMPFILE}"  "${OUTFILE3}"

/zhome/f5/1/58775/anaconda3/envs/mpa/bin/metaphlan "${BOWTIEFILE}" \
	-t marker_ab_table \
	--bowtie2db /work3/apca/metaphlan/databases \
    --nproc 8 \
    --input_type bowtie2out \
    > "${TEMPFILE}"

mv "${TEMPFILE}"  "${OUTFILE4}"

echo `date` Deleting the bowtie2 file

rm -f "${BOWTIEFILE}"

echo `date` Finished!