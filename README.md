## DSP_NF-METAGENOMICS
<p align="left">
Shotgun metagenomics pipeline to process microbiome samples
</p>

---

## Table of contents

- [About](#about) - Overview of the project's purpose and goals
- [Getting Started](#getting-started) - Instructions on how to begin with this project
- [Prerequisites and installing](#prerequisites-and-installing) - Required software and installation steps 
- [Step by step](#step-by-step) - Detailed guide to each stage of the project
- [Bioinformatic parameters](#bioinformatic-parameters) - Explanation and details of the bioinformatic parameters used throughout the pipeline
- [Repository structure](#repository-structure) - A layout of the repository's architecture, describing the purpose of each file or directory
- [References](#references) - Tools used in the project
- [Authors](#authors) - List of contributors to the project
- [Acknowledgments](#acknowledgement) - Credits and thanks to those who helped with the project

## About <a name = "about"></a>
The repository presents a comprehensive workflow for metagenomic analysis, starting from an initial assessment of data quality to an 
in-depth understanding of the composition and function of the examined microbiome. The analysis begins with a quality check of the 
sequenced data using FastQC, followed by a specific quality control for metagenomic data with Kneaddata. Subsequently, the workflow 
proceeds to the assembly of the reads with MegaHit and the classification of contigs into eukaryotic or prokaryotic. Anvi'o is then 
employed for the taxonomic and functional annotation of the contigs, as well as for mapping high-quality reads. Finally, Metaphlan 4.0 
facilitates further taxonomic annotation and the estimation of the abundance of various species based on reference genomes, thus 
completing the comprehensive analysis of the microbiome.

## Getting started <a name = "getting-started"></a>
The following instructions are designed to guide users in extracting information from their FASTQ files. Originally, the pipeline was implemented using 
shell scripts that invoke various bioinformatics software for data analysis. Presently, it is undergoing a transition to be re-implemented as a 
[Nextflow](https://nextflow.io) metagenomics workflow. This update aims to enhance the reproducibility and efficiency of the analysis process.

## Prerequisites and installing <a name = "prerequisites-and-installing"></a>
This workflow is configured to be executed through Docker and Azure Batch, leveraging cloud computing resources and containerized environments.
No additional installation steps are required for the workflow itself, as it relies on cloud-based resources and Docker containers.
Make sure Docker is installed and properly set up, configure your Azure Blob Storage and Azure Batch accounts, and install Nextflow following [Nextflow 
information guide](https://www.nextflow.io/docs/latest/getstarted.html) if you haven't done it yet. Once these prerequisites are in place, you can clone the 
repository and run the analysis.

Once you have cloned the repository from GitHub, it is important to configure the `nextflow.config` file, considering the paths of your container in your Azure account.
Then type `touch credentials.json` and copy and paste the private keys into the respective storage account and batch account names: `{
    "storageAccountName": "ma****ge",
    "storageAccountKey": "****Mi7MWBz****==",
    "batchAccountName": "***dtu***",
    "batchAccountKey": "****wX7rHYMD****=="
}`.
After that, enter the command `nextflow run main.nf -c nextflow.config -profile <name-profile-on-config-file> -w az://<your-container-name>` to run the pipeline.


## Step by step <a name = "step-by-step"></a>
1. Quick QC check of the raw sequenced data (fastQC)
2. quality control of metagenomic data, meant for microbiome experiments (Kneaddata)
3. Assembly of the reads (MegaHit): per sample and coassembly (step previous to Anvi'o)
4. Kingdom distribution: Prediction of whether a contig is eukaryotic or prokaryotic   
5. Taxonomical and functional annotation of the contigs (Anvi'o tools)
6. Mapping high-quality reads to the contigs (within Anvi'o framework)
7. Taxonomical annotation and taxa abundance estimation based on reference genomes (Metaphlan 4.0)

## Bioinformatic parameters <a name = "bioinformatic-parameters"></a>
Below is a detailed overview of the parameters used in each bioinformatic tool within the Nextflow pipeline (file: `nextflow_orange.nf`), specifically 
outlining the commands and their functions within the context of the entire workflow.

**FASTQC**
Tool designed for the quality control analysis og high-throughput sequencing data reporting visualizations that help assess the quality and characteristics of sequencing data before downstream analysis.
| Command         | Description                                                |
|-----------------| -----------------------------------------------------------| 
| -o (--output)   | Specifies the output directory to store the processed data.|
| -q              | Specifies the ... .                                        |

**KNEADDATA**
Tool used for QC and pre-processing of metagenomic and metatranscriptomic sequencing data;we need to consider we are working with input paired-end sequences files.
| Command         | Description                                                            |                                                                                                                                                 
|-----------------| -----------------------------------------------------------------------|                                        
| -i1             | Specifies the path to the input file containing the forward (R1) reads.|
| -i2             | Specifies the path to the input file containing the reverse (R2) reads.|
| --reference-db  | Specifies the reference database to be used for contaminant removal.   |
| --output        | Specifies the output directory to store the processed data.            |
| --bypass-trim   | Skip the trimming step during the processing of sequencing data.       |

**MEGAHIT**
Metagenome assembly tool used for assembling seqeuncing data particularly obtained from high-throughput sequencing technologies.
| Command         | Description                                                                         |
|-----------------| ------------------------------------------------------------------------------------|
| -1              | Specifies the path to the input file containing the first pair of paired-end reads. |
| -2              | Specifies the path to the input file containing the second pair of paired-end reads.|
| -o              | Specifies the output directory to store the assembled contigs or output files.      |

**WHOKARYOTE**
Tool which uses random forest to rpedict wheter a contig is from eukaryote or from a prokaryote(https://github.com/LottePronk/whokaryote).
| Command         | Description                                                                              |
|-----------------| -----------------------------------------------------------------------------------------|
| --contigs       | Specifies the path with your contigs file.                                               |
| --minsize       | Specifies a minimum contig size in bp, by default is 5000 (accuracy below 5000 is lower).|
| --outdir        | Specifies the output directory to store the output file.                                 |

**METAPHLAN**
Tool used for taxonomic profiling of metagenomic sequencing data (used for identification and quantification of microbial species present in a given sample based on unique clade-specific marker genes)  <br>
| Command         | Description                                                                                                                                         |
|-----------------| ----------------------------------------------------------------------------------------------------------------------------------------------------|
| -t              | Specifies the taxonomic level for the output; it allows users to choose the level of taxonomic resolution for the results.                          |
| --bowtie2out    | Specifies the output file for Bowtie2 alignments generated; it is used internally by MetaPhlAn for read alignments against the marker gene database.|
| --input_type    | Specifies the input data type for MetaPhlAn. It allows users to inform MetaPhlAn about the type of input data being provided (fastq, sam, bam).     |

## Repository structure <a name="repository-structure"></a>
The table below provides an overview of the key files and directories in this repository, along with a brief description of each.
|File  |Description            |
|:----:|-----------------------|
|[bin/](bin/)|Folder with python scripts adapted to the workflow|
|[map/](map/)|Folder with pdf and png for better rapresent the workflow|
|[old_scripts](old_scripts)|Folder with all the scripts used for creating the workflow (qc, assemblying, predictions, taxonimical annotation, mapping, etc...|
|[nextflow.config](nextflow.config)|Configuration file which contains a nextflow configuration for running the bioinformatics workflow, including parameters for processing genomic data on Azure cloud service|
|[nextflow_config_full_draft.txt](nextflow_config_full_draft.txt)|Text file which contains a configuration for nextflow workflow specifying resources requirements for each program used|

Add additional notes about how to deploy this on a live system.
## References <a name = references"></a>
- [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [KNEADDATA](https://huttenhower.sph.harvard.edu/kneaddata/)
- [MEGAHIT](https://www.metagenomics.wiki/tools/assembly/megahit)
- [WHOKARYOTE](https://github.com/LottePronk/whokaryote)
- [METAPHLAN](https://github.com/biobakery/MetaPhlAn)

## Authors <a name = "authors"></a
Contact me at apca@biosustain.dtu.dk if you are interested in running it before it is done.
- [apca@biosustain.dtu.dk](https://github.com/apalleja)
- [marcor@dtu.dk](https://github.com/marcoreverenna)

## Acknowledgements <a name = "acknowledgement"></a>
We would like to extend our heartfelt gratitude to [DTU Biosustain](https://www.biosustain.dtu.dk/) and the Novo Nordisk Foundation 
Center for Biosustainability for providing the essential resources and support that have been 
fundamental in the development and success of the [DSP (Data Science 
Platform)](https://www.biosustain.dtu.dk/informatics/data-science-platform) and [MoNA (Multi-omics Network 
Analysis)](https://www.biosustain.dtu.dk/research/research-groups/multi-omics-network-analytics-alberto-santos-delgado) projects.


