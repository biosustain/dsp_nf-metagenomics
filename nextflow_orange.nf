/*
 * Pipeline input parameters
 */
params.reads = "$projectDir/data/O*_{1,2}.fq.gz"
params.orange_genome = "$projectDir/orange_genome/GCF_022201045.2_DVS_A1.0_genomic.fna"
params.multiqc = "$projectDir/multiqc"
params.genomedir = "orange_genome"
params.outdir = "results"
log.info """\
    ORANGE PEEL METAGENOMICS - N F   P I P E L I N E
    ===================================
    orange_genome: ${params.orange_genome}
    reads        : ${params.reads}
    genomedir    : ${params.genomedir}      
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
    publishDir params.outdir, mode: "copy"

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

    tag "MultiQC on all fastQC"
    //publishDir params.outdir, mode:'copy'

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

    tag "Kneaddata on $sample_id"

    input:
    path orange_genome
    tuple val(sample_id), path(reads)

    output:
    //path "$sample_id"
    path("*_1_kneaddata_paired_{1,2}.fastq")
    //, emit: kneaddata_qc

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
    container "quay.io/biocontainers/megahit:1.2.9--h43eeafb_4"

    tag "Megahit"
    publishDir "./assembly", mode:'copy'
    cpus 4

    input:
    path(files)
    //set sample_id, file(files) from qc_ch

    output:
    path("*.contigs.fa")

    script:
    """
    megahit -1 ${files[0]} -2 ${files[1]} \
	-t $task.cpus \
    -o output
    """

}

/*
 * Whokaryote to predict wheter the contigs are eukaryotic or prokaryotic
 */
process WHOKARYOTE {
    container "quay.io/biocontainers/whokaryote"

    tag "Whokaryote on contigs of $sample_id"

    input:
    tuple val(sample_id), path(files)

    output:
    path "${sample_id}.contigs.fa"

    script:
    """
    whokaryote --contigs ${files[0]} \
	--minsize 1000 \
	--outdir .
    """
}

/*
 * Running Metaphlan on the high quality reads
 */
process METAPHLAN {
    container "quay.io/biocontainers/metaphlan"

    tag "Metaphlan on HQ reads of $sample_id"
    cpus 4

    input:
    tuple val(sample_id), path(files)

    output:
    path "metaphlan_abundance.tsv"

    script:
    """
    metaphlan "${files[0]},${files[1]}" \
	-t rel_ab \
    --bowtie2out . \
    --nproc $task.cpus \
    --input_type fastq \
    > .
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

    //read_pairs_ch.view()
    //return
    index_ch = BUILD_HOST_DB(params.orange_genome)
    //index_ch.view()
    //return

    fastqc_ch = FASTQC(read_pairs_ch)
    //index_ch.view()
    qc_ch = QC(index_ch, read_pairs_ch)
    MULTIQC(qc_ch.mix(fastqc_ch).collect())
    qc_ch.view()
    assembly_ch = ASSEMBLY(qc_ch)
    assembly_ch.view()
    return

    who_ch = WHOKARYOTE(assembly_ch)
    mpa_ch = METAPHLAN(qc_ch)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
