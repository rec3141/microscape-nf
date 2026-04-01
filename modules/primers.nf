// Primer removal with cutadapt.
//
// Removes primer sequences from paired-end reads using separate forward and
// reverse primer files. R1 must match a forward primer and R2 must match a
// reverse primer — pairs missing either are discarded.
//
// Auto-detection selects primer set by filename prefix (16S, 18S, ITS).
// Default: bacterial 16S (515F/806RB).

process REMOVE_PRIMERS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/trimmed", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/trimmed" : null

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1.trimmed.fastq.gz"), path("${meta.id}_R2.trimmed.fastq.gz"), emit: reads
    path("${meta.id}_cutadapt.log"), emit: log

    script:
    // Auto-select primer files based on sample name prefix
    def fwd_file = params.primers_fwd ?: (
        params.primer_auto ? (
            meta.id =~ /^18/ ? params.primers_fwd_euk :
            meta.id =~ /^ITS/ ? params.primers_fwd_its :
            params.primers_fwd_bac
        ) : params.primers_fwd_bac
    )
    def rev_file = params.primers_rev ?: (
        params.primer_auto ? (
            meta.id =~ /^18/ ? params.primers_rev_euk :
            meta.id =~ /^ITS/ ? params.primers_rev_its :
            params.primers_rev_bac
        ) : params.primers_rev_bac
    )
    """
    cutadapt \\
        -g file:${fwd_file} \\
        -G file:${rev_file} \\
        --discard-untrimmed \\
        --pair-filter=any \\
        -j ${task.cpus} \\
        -e ${params.primer_error_rate} \\
        -o ${meta.id}_R1.trimmed.fastq.gz \\
        -p ${meta.id}_R2.trimmed.fastq.gz \\
        ${r1} ${r2} \\
        > ${meta.id}_cutadapt.log 2>&1

    echo "[INFO] ${meta.id}: primer removal complete" >&2
    """
}
