#! /bin/bash
#
#  create_sampletable.sh
#  2023-08-08
#
#  Finds all sequencing data within the given path and creates a "sampletable",
#  in which duplicated filenames are represented only once.  Fastq files must be to
#  be paired end, with filenames ending in "_1.fq.gz" and "_2.fq.gz".
#
#  Sample IDs (first column of the sampletable) are inferred from the parent directory
#  name, because this is the file structure currently used by for example Novogene.
#
#  This script probably works with data provided by different sequencing providers, BUT CHEK IT OUT!.
#
#  The recommended usage below creates a log file "sampletable.log".
#
#  Usage:
#
#    bash create_novogene_sampletable.sh  <path>  |& tee sampletable.log
#
#  A concrete example:
#
#    bash create_novogene_sampletable.sh  /home/projects/cu_10051/data/nesifo  |& tee sampletable.log
#

# Stop with an error if anything fails
set -euo pipefail

TEMP1="create_novogene_sampletable_1.TEMP"
TEMP2="create_novogene_sampletable_2.TEMP"
OUTFILE="sampletable"

echo "Creating sampletable $(date -Is)"
echo "Searching for files in $*"

# First find all forward read files
find $@ -name '*_1.fq.gz' -print > ${TEMP1}

echo "Found $(cat ${TEMP1} | wc -l) files ending in _1.fq.gz"

# Now remove duplicated filenames (basenames), even if they are in different directories
paste ${TEMP1}  <(cat ${TEMP1} | xargs -n 1 basename)  \
| sort -k 2 | uniq -f 1 | cut -f 1  \
> ${TEMP2}

echo "  of which $(cat ${TEMP2} | wc -l) files were not duplicates and thus retained"

# Create full sampletable (samplename, fq1, fq2) based on these unique filenames
paste <(perl -lpe 's|^.+/([^/]+)/[^/]+$|\1|' ${TEMP2})  \
      ${TEMP2}  \
      <(sed 's/_1.fq.gz/_2.fq.gz/' ${TEMP2})  \
> ${OUTFILE}

# Clean up temporary files
rm ${TEMP1} ${TEMP2}

echo "Finished with sampletable!"