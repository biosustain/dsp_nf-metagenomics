## DSP_NF-METAGENOMICS
<p align="left">
Shotgun metagenomics pipeline to process microbiome samples
</p>

---

## 📝 Table of contents

- [About](#about)
- [Getting Started](#getting_started)
- [Prerequisites and installing](#prerequisites_and_installing)
- [Deployment](#deployment)
- [Usage](#usage)
- [Built Using](#built_using)
- [Authors](#authors)
- [Acknowledgments](#acknowledgement)

## 🧐 About <a name = "about"></a>
The objective of this project is ...

## 🏁 Getting Started <a name = "getting_started"></a>
These instructions will enable you to ...
This is a pipeline that was written initially in shell scripts that call different software and in few python and R scripts.
It is currently in the process to be coded as a nextflow metagenomics workflow.

## 🔧 Prerequisites and installing <a name = "prerequisites_and_installing"></a>
This workflow is set up to be executed on ...
Therefore, the only prerequisite is to have ...

1. Install dependencies into isolated environment
```
conda env create --name your_conda_name --file file.yaml
```
2. Activate environment
```
source activate imp_proj
```

## 🚀 Deployment <a name = "deployment"></a>
1. Quick QC check of the raw sequenced data (fastQC)
2. quality control of metagenomic data, meant for microbiome experiments (Kneaddata)
3. Assembly of the reads (MegaHit): per sample and coassembly (step previous to Anvi'o)
4. Kingdom distribution: Prediction of whether a contig is eukaryotic or prokaryotic   
5. Taxonomical and functional annotation of the contigs (Anvi'o tools)
6. Mapping high-quality reads to the contigs (within Anvi'o framework)
7. Taxonomical annotation and taxa abundance estimation based on reference genomes (Metaphlan 4.0)

## 🧬 Bioinformatic parameters <a name = "deployment"></a>
Following, the parameters included in the nextflow pipeline (file --> nextflow_orange.nf)

1. FASTQC
short description: tool designed for the quality control analysis og high-throughput sequencing data reporting visualizations that help assess the quality and 
characteristics of sequencing data before downstream analysis
a. -o (--output) = specifies the output directory to store the processed data
b. -q = ?

2. KNEADDATA
short description: tool used for QC and pre-processing of metagenomic and metatranscriptomic sequencing data;we need to consider we are working with input paired-end 
sequences files.
a. -i1 = specifies the path to the input file containing the forward (R1) reads.
b. -i2 = specifies the path to the input file containing the reverse (R2) reads.
c. --reference-db = specifies the reference database or databases to be used for contaminant removal (host genomes, known contaminants, etc...)
d. --output = specifies the output directory to store the processed data
e. --bypass-trim = skip the trimming step during the processing of sequencing data (either already trimmed or not required)

3. MEGAHIT
short description: metagenome assembly tool used for assembling seqeuncing data particularly obtained from high-throughput sequencing technologies
a. -1 = specifies the path to the input file containing the first pair of paired-end reads.
b. -2 = specifies the path to the input file containing the second pair of paired-end reads.
c. -o = specifies the output directory to store the assembled contigs or output files

4. WHOKARYOTE
short description: https://github.com/LottePronk/whokaryote; this tool uses random forest to rpedict wheter a contig is from eukaryote or from a prokaryote
a. --contigs = specifies the path with your contigs file
b. --minsize = specifies a minimum contig size in bp, by default is 5000 (accuracy on contigs below 5000 is lower)
c. --outdir = specifies the output directory to store the output file

5. METAPHLAN
short description: tool used for taxonomic profiling of metagenomic sequencing data (used for identification and quantification of microbial species present in a 
given sample based on unique clade-specific marker genes)
a. -t = specifies the taxonomic level for the output; it allows users to choose the level of taxonomic resolution for the results
b. --bowtie2out = specifies the output file for Bowtie2 alignments generated --> it is used internally by MetaPhlAn for read alignments against the marker gene database. 
This parameter allows users to specify the file path where the Bowtie2 alignment results will be stored.
c. --input_type = specifies the input data type for MetaPhlAn. It allows users to inform MetaPhlAn about the type of input data being provided (fastq, sam, bam) 


## 🎈 Usage <a name="usage"></a>
Add notes about how to use the system.

Add additional notes about how to deploy this on a live system.
## ⛏️ Built Using <a name = "built_using"></a>
- [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [KNEADDATA](https://huttenhower.sph.harvard.edu/kneaddata/)
- [MEGAHIT](https://www.metagenomics.wiki/tools/assembly/megahit)
- [WHOKARYOTE](https://github.com/LottePronk/whokaryote)
- [METAPHLAN](https://github.com/biobakery/MetaPhlAn)
## ✍️ Authors <a name = "authors"></a>
- [apca@biosustain.dtu.dk](https://github.com/apalleja)
Contact me at apca@biosustain.dtu.dk if you are interested in running it before it is done.
## 🎉 Acknowledgements <a name = "acknowledgement"></a>
- MoNA
- DTU
- Novo Nordisk Foundation

