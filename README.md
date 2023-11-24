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

## 🎈 Usage <a name="usage"></a>
Add notes about how to use the system.

Add additional notes about how to deploy this on a live system.
## ⛏️ Built Using <a name = "built_using"></a>
- [GenomeDK](https://genome.au.dk/) - Server Environment
- [Beagle5.4](https://faculty.washington.edu/browning/beagle/beagle.html) - Imputation algorithm
## ✍️ Authors <a name = "authors"></a>
- [apca@biosustain.dtu.dk](https://github.com/apalleja)
Contact me at apca@biosustain.dtu.dk if you are interested in running it before it is done.
## 🎉 Acknowledgements <a name = "acknowledgement"></a>
Add notes...
