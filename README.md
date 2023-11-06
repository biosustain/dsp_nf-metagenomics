# dsp_nf-metagenomics
Shotgun metagenomics pipeline to process microbiome samples

It has different steps:

1. Quick QC check of the raw sequenced data (fastQC)
2. quality control of metagenomic data, meant for microbiome experiments (Kneaddata)
3. Assembly of the reads (MegaHit): per sample and coassembly (step previous to Anvi'o)
4. Kingdom distribution: Prediction of whether a contig is eukaryotic or prokaryotic
5. Taxonomical and functional annotation of the contigs (Anvi'o tools)
6. Mapping high-quality reads to the contigs (within Anvi'o framework)
7. Taxonomical annotation and taxa abundance estimation based on reference genomes (Metaphlan 4.0)

This is a pipeline that was written initially in shell scripts that call different software and in few python and R scripts. It is currently in the process to be coded as a nextflow metagenomics workflow. Contact me at apca@biosustain.dtu.dk if you are interested in running it before it is done.
