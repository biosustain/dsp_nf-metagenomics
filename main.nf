/*
 * Pipeline input parameters
 */
//params.reads = "$projectDir/data/*_{1,2}.fq.gz"
//params.host_genome = "$projectDir/host_genome/GCF_022201045.2_DVS_A1.0_genomic.fna"
//params.genomedir = "$projectDir/host_genome"
//params.metaphlan_db = "$projectDir/databases/metaphlan_db"
//params.eggnog_db = "$projectDir/databases/eggnog_db"
//params.outdir = "$params.outdir/results2"

println "projectDir: $projectDir"

log.info """\
    METAGENOMICS - N F   P I P E L I N E
    ===================================
    host_genome  : ${params.host_genome}
    reads        : ${params.reads}
    genomedir    : ${params.genomedir}
    metaphlan_db : ${params.metaphlan_db}
    eggnog_db    : ${params.eggnog_db}
    outdir       : ${params.outdir}
    """
    .stripIndent(true)

/*
 * Preprocess fastq files
 */
process PREPROCESS_FASTQ {
    container "quay.io/bedrock/ubuntu"
    publishDir "${params.outdir}/data_processed", mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_{1,2}.fq.gz"), emit: processed_data

    script:
    """
    gunzip -c ${reads[0]} | sed 's/\\(^@[^ ]*\\) .*/\\1\\/1/' | gzip > ${sample_id}_1.fq.gz
    gunzip -c ${reads[1]} | sed 's/\\(^@[^ ]*\\) .*/\\1\\/2/' | gzip > ${sample_id}_2.fq.gz
    """
}

/*
 * Building a database for the orange genome
 */
process BUILD_HOST_DB {
    container "quay.io/biocontainers/bowtie2:2.5.1--py36hb79b6da_1"
    publishDir params.genomedir, mode: "copy"
    cpus 4

    input:
    path host_genome

    output:
    path "host_genome.*"

    script:
    """
    bowtie2-build --threads $task.cpus ${host_genome} host_genome
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
 * Using Adapter removal to QC the raw reads
*/
process ADAPTERREMOVAL {
    tag "Adapter Removal - $sample_id"

    container "quay.io/biocontainers/adapterremoval:2.3.2--hb7ba0dd_0"
    cpus = 10

    publishDir "${params.outputDir}/adapterremoval_logs/", pattern: "${sample_id}.{log,progress}", mode:'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.fq.gz"), emit: trimmed_reads
    path "${sample_id}.log", emit: adapter_log
    path "${sample_id}.progress.log"

    script:
    """
    AdapterRemoval                                                      \
        --threads ${task.cpus}                                          \
        ${params.adapterTrimParameters}                                 \
        --file1 ${reads[0]}                                             \
        --file2 ${reads[1]}                                             \
        --output1 "${sample_id}_hf_trimmed_1.fq.gz"                     \
        --output2 "${sample_id}_hf_trimmed_2.fq.gz"                     \
        --gzip  --gzip-level 1                                          \
        --settings "${sample_id}.log"                                   \                                          \
        > >(tee -a "${sample_id}.progress.log")                         \
        2> >(tee -a "${sample_id}.progress.log" >&2)
    """
}

/*
 * Host removal
 */
process HOSTFILTER {
    tag "Hostfilter - $sample_id"

    container "quay.io/biocontainers/bowtie2:2.5.1--py36hb79b6da_1"
    cpus = 16

    publishDir "${params.outputDir}/hostfilter_bowtie_logs/", pattern: "log/bowtie2/${sample_id}.log", mode:'copy'
    publishDir "${params.outputDir}/hostfilter_samtools_logs/", pattern: "log/samtools/${sample_id}.log", mode:'copy'

    input:
    path index_files
    tuple val(sample_id), path(reads_1), path(reads_2)

    output:
    tuple val(sample_id), path("*.fq.gz"),  emit: filtered_reads
    path "hostfilter_bowtie_logs/${sample_id}.log", emit: bowtie2_log
    path "hostfilter_samtools_logs/${sample_id}.log", emit: samtools_log

    script:
    def idx = file(params.reference[ params.geneCatalogName ].hostDatabase).getName()
    """
    mkdir -p log/bowtie2 log/samtools
    bowtie2                                                     \
            -x ${idx}                                           \
            -1 ${reads_1.join(',')}                             \
            -2 ${reads_2.join(',')}                             \
            -k 1                                                \
            --threads ${task.cpus}                              \
            --phred33                                           \
            --local                                             \
            2> >(tee -a "log/bowtie2/${sample_id}.log" >&2)     \
        | samtools fastq                                        \
            -f 12  -F 2304  -n                                  \
            -1 >( pigz -p 6 > ${sample_id}.hf_1.fq.gz )         \
            -2 >( pigz -p 6 > ${sample_id}.hf_2.fq.gz )         \
            -0 /dev/null                                        \
            -s /dev/null                                        \
            2> >(tee -a "log/samtools/${sample_id}.log" >&2)
    sleep 2   # make sure pigz processes finish
    """
}
/*
 * Quality control of the reads, decontamination
 */
process QC {
    
    container "biobakery/kneaddata"
    tag "Kneaddata on $sample_id"
    cpus 8

    input:
    path host_genome
    tuple val(sample_id), path(processed_data)

    output:
    tuple val(sample_id), path("${sample_id}_1_kneaddata_paired_{1,2}.fastq"), emit: kneaddata_qc
    path "kneaddata_logs/${sample_id}_1_kneaddata.log", emit: kneaddata_logs
    path "kneaddata_fastqc_end/${sample_id}_EKDN24000116{3,4}-1A_222VF7LT4_L{4,5}_fastqc.{html,zip}", emit: fastqc_end
    
    script:
    """
    kneaddata --input ${processed_data[0]} --input ${processed_data[1]} \
    -t $task.cpus \
	-db host_genome \
	--output . \
    --run-fastqc-end \
    --trimmomatic-options "$params.trimmomatic_options"

    mkdir -p kneaddata_logs
    mv ${sample_id}_1_kneaddata.log kneaddata_logs/
    """
}

/*
 * Create a table of QC counts for all samples
 */
process QC_STATS {
    container "quay.io/biocontainers/kneaddata:0.12.0--pyhdfd78af_1"
    tag "Kneaddata read count for all samples"

    input:
    path("kneaddata_logs/*_1_kneaddata.log") // it worked with *.log, also like that but problem with sample name

    output:
    path "kneaddata_read_count_table.tsv"
    
    script:
    """
    kneaddata_read_count_table --input kneaddata_logs \
    --output kneaddata_read_count_table.tsv
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

    fastqc_ch = FASTQC(read_pairs_ch)
    fastqc_ch.view()
    MULTIQC(fastqc_ch.collect())
    return

    index_ch = BUILD_HOST_DB(params.host_genome)
    index_ch.view()

    
    

    index_ch = BUILD_HOST_DB(params.host_genome)
    index_ch.view()
    
    QC(index_ch, fastq_processed_ch)
    QC.out.kneaddata_qc.view()
    //logFilesChannel = Channel.Create()
    // { file -> logFilesChannel.put(file) }
    QC_STATS(QC.out.kneaddata_logs.collect())

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
