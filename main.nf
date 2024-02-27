/*
 * Pipeline input parameters
 */
params.reads = "$projectDir/data/*_{1,2}.fq.gz"
params.orange_genome = "$projectDir/orange_genome/GCF_022201045.2_DVS_A1.0_genomic.fna"
params.multiqc = "$projectDir/multiqc"
params.genomedir = "$projectDir/orange_genome"
params.metaphlan_db = "$projectDir/databases/metaphlan_db"
params.eggnog_db = "$projectDir/databases/eggnog_db"
params.outdir = "$params.outdir/results2"

println "projectDir: $projectDir"

log.info """\
    METAGENOMICS - N F   P I P E L I N E
    ===================================
    orange_genome: ${params.orange_genome}
    reads        : ${params.reads}
    genomedir    : ${params.genomedir}
    metaphlan_db : ${params.metaphlan_db}
    eggnog_db    : ${params.eggnog_db}
    outdir       : ${params.outdir}
    """
    .stripIndent(true)

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
    publishDir params.outdir, mode:'copy'

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
    tag "Kneaddata on $sample_id"

    input:
    path orange_genome
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_1_kneaddata_paired_{1,2}.fastq"), emit: kneaddata_qc
    tuple val(sample_id), path("${sample_id}_1_kneaddata.log"), emit: kneaddata_log
    
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
    publishDir "${params.outdir}/assembly", mode:'copy'
    tag "Megahit on $sample_id"
    
    input:
    tuple val(sample_id), path(kneaddata_qc)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}.contigs.fa"), emit: contigs_id
    path("${sample_id}/${sample_id}.contigs.fa"), emit: contigs
    path("${sample_id}/${sample_id}.log"), emit: contigs_logs

    script:
    """
    megahit -1 ${kneaddata_qc[0]} -2 ${kneaddata_qc[1]} \
    --out-prefix ${sample_id} \
    -o ${sample_id}
    """
}

/*
 * Summarising Assembly results
 */
process ASSEMBLY_STATS {
    container 'pandas/pandas:pip-all'
    publishDir "${params.outdir}/assembly", mode:'copy'
    tag "Summarising Assembly results"

    input:
    path(contigs)
    path(contigs_logs)

    output:
    path "contigs_summary.txt"
    path "assembly_stats.txt"

    script:
    """
    dist_contig_lengths.py -a $contigs -o contigs_summary.txt
    assembly_stats.py -a $contigs_logs -o assembly_stats.txt
    """
}

/*
 * Calling genes on the contigs with Prodigal
 */
 process CALLGENES {
    container "quay.io/biocontainers/prodigal:2.6.3--h516909a_2"
    publishDir "${params.outdir}/prodigalGenes", mode:'copy'
    tag "Calling genes on $sample_id contigs"
    
    input:
    tuple val(sample_id), path(contigs_id)
    
    output:
    tuple val(sample_id), path("${sample_id}.gff"), emit: gene_annotations
    tuple val(sample_id), path("${sample_id}.fna"), emit: nucleotide_fasta
    tuple val(sample_id), path("${sample_id}.faa"), emit: amino_acid_fasta
    tuple val(sample_id), path("${sample_id}_all.txt"), emit: all_gene_annotations

    script:
    """
    prodigal -f gff -d ${sample_id}.fna \
    -o ${sample_id}.gff \
    -a ${sample_id}.faa \
    -s ${sample_id}_all.txt \
    -i ${contigs_id}
    """
}

/*
 * Whokaryote to predict the contigs kingdom
 */
process WHOKARYOTE {
    container "quay.io/biocontainers/whokaryote:1.0.1--pyhdfd78af_0"
    publishDir "${params.outdir}/whokaryote", mode:'copy'
    tag "Whokaryote on contigs of $sample_id"

    input:
    tuple val(sample_id), path(contigs_id)

    output:
    path("${sample_id}/${sample_id}.featuretable_predictions_T.tsv"), emit: predictions

    script:
    """
    whokaryote.py --contigs ${contigs_id} \
	--minsize 1000 \
	--outdir ${sample_id}

    mv "${sample_id}/featuretable_predictions_T.tsv" "${sample_id}/${sample_id}.featuretable_predictions_T.tsv"
    """
}

/*
 * Summarising Whokaryote results
 */
process WHOKARYOTE_STATS {
    container 'pandas/pandas:pip-all'
    publishDir "${params.outdir}/whokaryote", mode:'copy'
    tag "Summarising Whokaryote results"

    input:
    path(predictions)

    output:
    path "whokaryote_stats.txt"

    script:
    """
    whokaryote_stats.py -p $predictions -o whokaryote_stats.txt
    """
}

/*
 * Running Metaphlan on the high quality reads
 */
process METAPHLAN {
    container "quay.io/biocontainers/metaphlan:4.0.6--pyhca03a8a_0"
    publishDir "${params.outdir}/metaphlan", mode:'copy'
    cpus 8

    tag "Metaphlan on HQ reads of $sample_id"

    input:
    tuple val(sample_id), path(kneaddata_qc)

    output:
    path("${sample_id}.profiled_metagenome.txt")
    

    script:
    """
    rm -rf bowtie
    mkdir bowtie

    metaphlan "${kneaddata_qc[0]},${kneaddata_qc[1]}" \
	-t rel_ab \
    --bowtie2db ${params.metaphlan_db} \
    --bowtie2out bowtie/${sample_id}.bowtie2.bz2 \
    --nproc $task.cpus \
    --input_type fastq \
    --output_file ${sample_id}.profiled_metagenome.txt
    """  
}

/*
 * Merging metaphlan abundance tables
 */
 process METAPHLAN_MERGE {
    container "quay.io/biocontainers/metaphlan:4.0.6--pyhca03a8a_0"
    publishDir "${params.outdir}/metaphlan", mode:'copy'

    input:
    path(files)

    output:
    path "merged_abundance_table.txt"

    script:
    """
    merge_metaphlan_tables.py ${files} \
    -o merged_abundance_table.txt
    """  
 } 

/*
 * Downloading EggNOG database
 */
 process EGGNOG_DB_DOWNLOAD {

    container "metashot/eggnog-mapper:2.1.4-2"
    publishDir "${params.outdir}/databases" , mode: 'copy'
    tag 'Downloading EggNOG database'

    output:
    path 'eggnog_db', type: 'dir', emit: eggnog_db

    script:
    """
    mkdir eggnog_db
    download_eggnog_data.py -y --data_dir eggnog_db
    """
}


/*
 * Mapping read with EggNOG
 */
 process EGGNOG_MAPPER {
    //container "quay.io/biocontainers/eggnog_mapper:2.1.12--pyhdfd78af_0"
    container "metashot/eggnog-mapper:2.1.4-2"
    publishDir "${params.outdir}/eggnog" , mode: 'copy'
    tag "EGGNOG mapping on $sample_id"
    cpus 8

    input:
    tuple val(sample_id), path(proteins)
    path(eggnog_db)
   
    output:
    path "${sample_id}*"
    path "${sample_id}.emapper.annotations", emit: annotations
   
    script:
    param_eggnog_db_mem = params.eggnog_db_mem ? '--dbmem' : ''
    """
    mkdir temp

    emapper.py \
        -i ${proteins} \
        -o ${sample_id} \
        -m diamond \
        --itype proteins \
        --temp_dir temp \
        --data_dir ${eggnog_db} \
        --cpu ${task.cpus} \
        ${param_eggnog_db_mem}
        
    rm -rf temp
    """
}

/*
 * Merging EGGNOG mapper annotations
 */
 process MERGE_EGGNOG_MAPPER {
    container = 'metashot/utils:1.1.0-2'
    publishDir "${params.outdir}/eggnog_tables" , mode: 'copy'
    tag "Merging EGGNOG annotations"

    input:
    path(annotations)
   
    output:
    path 'eggnog_*.tsv'

    script:
    """
    merge_eggnog_mapper.py -i ${annotations}
    """
}

/*
 * Definning the workflow
 */
workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .toSortedList( { a, b -> a[0] <=> b[0] } )
        .flatMap()
        .set { read_pairs_ch }
        read_pairs_ch.view()

    index_ch = BUILD_HOST_DB(params.orange_genome)
    index_ch.view()

    fastqc_ch = FASTQC(read_pairs_ch)
    fastqc_ch.view()
    MULTIQC(fastqc_ch.collect())

    QC(index_ch, read_pairs_ch)
    QC.out.kneaddata_qc.view()
    QC.out.kneaddata_log.view()

    ASSEMBLY(QC.out.kneaddata_qc)
    ASSEMBLY.out.contigs_id.view()
    ASSEMBLY.out.contigs.collect().view()
    ASSEMBLY.out.contigs_logs.collect().view()
    ASSEMBLY_STATS(ASSEMBLY.out.contigs.collect(), ASSEMBLY.out.contigs_logs.collect())
    
    who_ch = WHOKARYOTE(ASSEMBLY.out.contigs_id)
    who_ch.predictions.view()
    WHOKARYOTE_STATS(WHOKARYOTE.out.predictions.collect())

    mpa_cph = METAPHLAN(QC.out.kneaddata_qc)
    mpa_cph.view()
    METAPHLAN_MERGE(mpa_cph.collect())
    
    genescalled_ch = CALLGENES(ASSEMBLY.out.contigs_id)
    CALLGENES.out.gene_annotations.view()
    CALLGENES.out.nucleotide_fasta.view()
    CALLGENES.out.amino_acid_fasta.view()
    CALLGENES.out.all_gene_annotations.view()

    EGGNOG_DB_DOWNLOAD()
    EGGNOG_DB_DOWNLOAD.out.eggnog_db.view()
    EGGNOG_MAPPER(CALLGENES.out.amino_acid_fasta, EGGNOG_DB_DOWNLOAD.out.eggnog_db)
    EGGNOG_MAPPER.out.annotations.view()
    MERGE_EGGNOG_MAPPER(EGGNOG_MAPPER.out.annotations.collect())
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
