/*
 * Pipeline input parameters
 */
params.reads = "$projectDir/data/O*_{1,2}.fq.gz"
params.orange_genome = "$projectDir/orange_genome/GCF_022201045.2_DVS_A1.0_genomic.fna"
params.multiqc = "$projectDir/multiqc"
params.genomedir = "$projectDir/orange_genome"
params.metaphlan_db = "$projectDir/databases/metaphlan_db"
params.outdir = "$params.outdir/results"

log.info """\
    ORANGE PEEL METAGENOMICS - N F   P I P E L I N E
    ===================================
    orange_genome: ${params.orange_genome}
    reads        : ${params.reads}
    genomedir    : ${params.genomedir}
    metaphlan_db : ${params.metaphlan_db}      
    outdir       : ${params.outdir}
    """
    .stripIndent()

/*
 * Building a database for the orange genome
 */
process BUILD_HOST_DB {
    container "quay.io/biocontainers/bowtie2:2.5.1--py36hb79b6da_1"
    publishDir params.genomedir, mode: "copy"
    cpus 4

    input:
    path orange_genome

    output:
    path "orange_genome.*"

    script:
    """
    bowtie2-build --threads $task.cpus ${orange_genome} orange_genome
    """
}

/*
 * Quickly checking raw reads quality
 */
process FASTQC {
    container "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
    //publishDir params.outdir, mode: "copy"
    //label "process_medium"

    tag "FASTQC on $sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "fastqc_${sample_id}_logs"

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -q ${reads}
    """
}

/*
 * Preparing a multi QC report with all fastQC reports
 */
process MULTIQC {
    container "quay.io/biocontainers/multiqc:1.16--pyhdfd78af_0"
    //label "process_single"

    tag "MultiQC on all fastQC"

    input:
    path '*'

    output:
    path 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}

/*
 * Quality control of the reads, decontamination
 */
process QC {
    container "quay.io/biocontainers/kneaddata:0.12.0--pyhdfd78af_1"
    //label "process_high"
    tag "Kneaddata on $sample_id"

    input:
    path orange_genome
    tuple val(sample_id), path(reads)

    output:
    //path("${sample_id}_1_kneaddata_paired_{1,2}.fastq"), emit: kneaddata_qc
    path("*_1_kneaddata_paired_{1,2}.fastq")
    
    script:
    """
    kneaddata -i1 ${reads[0]} -i2 ${reads[1]} \
	--reference-db orange_genome \
	--output . \
	--bypass-trim
    """
}

/*
 * Assembly of the high quality reads
 */
process ASSEMBLY {
    container "quay.io/biocontainers/megahit:1.2.9--h5b5514e_3"
    //label "process_high"

    tag "Megahit"
    publishDir params.outdir, mode:'copy'

    input:
    path(files)
    tuple val(sample_id), path(reads)

    output:
    path("${sample_id}/final.contigs.fa")

    script:
    """
    megahit -1 ${files[0]} -2 ${files[1]} \
    -o ${sample_id}
    """

}

/*
 * Whokaryote to predict wheter the contigs are eukaryotic or prokaryotic
 */
process WHOKARYOTE {
    container "quay.io/biocontainers/whokaryote:1.0.1--pyhdfd78af_0"
    publishDir "${params.outdir}/whokaryote", mode:'copy'
    //label "process_single"
    tag "Whokaryote on contigs of $sample_id"

    input:
    path(files)
    tuple val(sample_id), path(reads)

    output:
    path "whokaryote/${sample_id}/featuretable_predictions_T.tsv"

    script:
    """
    mkdir whokaryote
    whokaryote.py --contigs ${files[0]} \
	--minsize 1000 \
	--outdir whokaryote/${sample_id}
    """
}


/*
 * Running Metaphlan on the high quality reads
 */
process METAPHLAN {
    container "quay.io/biocontainers/metaphlan:4.0.6--pyhca03a8a_0"
    publishDir "${params.outdir}/metaphlan", mode:'copy'
    cpus 4
    //label "process_high"
    //mv metaphlan/${sample_id}.incomplete.txt metaphlan/${sample_id}.profiled_metagenome.tsv
    //> metaphlan/${sample_id}.incomplete.txt

    tag "Metaphlan on HQ reads of $sample_id"

    input:
    path(files)
    tuple val(sample_id), path(reads)

    output:
    path "metaphlan/${sample_id}.profiled_metagenome.tsv"

    script:
    """
    rm -rf metaphlan bowtie
    mkdir metaphlan
    mkdir bowtie

    metaphlan "${files[0]},${files[1]}" \
	-t rel_ab \
    --bowtie2db ${params.metaphlan_db} \
    --bowtie2out bowtie/${sample_id}.bowtie2.bz2 \
    --nproc $task.cpus \
    --input_type fastq \
    --output_file metaphlan/${sample_id}.profiled_metagenome.tsv
    """  
}


/*
 * Definning the workflow
 */
workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs_ch }
        read_pairs_ch.view()

    index_ch = BUILD_HOST_DB(params.orange_genome)
    index_ch.view()

    fastqc_ch = FASTQC(read_pairs_ch)
    fastqc_ch.view()

    qc_ch = QC(index_ch, read_pairs_ch)
    MULTIQC(qc_ch.mix(fastqc_ch).collect())
    qc_ch.view()

    assembly_ch = ASSEMBLY(qc_ch, read_pairs_ch)
    assembly_ch.view()

    who_ch = WHOKARYOTE(assembly_ch, read_pairs_ch)
    mpa_ch = METAPHLAN(qc_ch, read_pairs_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
