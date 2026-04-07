// Primer removal with cutadapt.
//
// Two-pass approach:
//   1. DETECT_PRIMERS: runs all primer pairs on each sample, picks the best match
//   2. REMOVE_PRIMERS: runs the selected primer pair with --discard-untrimmed
//
// If metadata provides a primer_pair column, DETECT_PRIMERS is skipped.

process DETECT_PRIMERS {
    tag "${meta.id}"
    cpus 1
    conda "${projectDir}/envs/python.yml"

    input:
    tuple val(meta), path(r1), path(r2)
    path(primer_files)

    output:
    tuple val(meta), path(r1), path(r2), env(BEST_PRIMER), emit: detected

    script:
    """
    # Run cutadapt with each primer file, count passing reads
    BEST_PRIMER=""
    BEST_COUNT=0

    for pf in ${primer_files}; do
        COUNT=\$(cutadapt -g file:\$pf -G file:\$pf --discard-untrimmed \
            -j ${task.cpus} -e ${params.primer_error_rate} \
            -o /dev/null -p /dev/null \
            ${r1} ${r2} 2>&1 | grep "Pairs written" | grep -oP '[\\d,]+' | head -1 | tr -d ',')
        COUNT=\${COUNT:-0}
        echo "\$pf: \$COUNT pairs"
        if [ "\$COUNT" -gt "\$BEST_COUNT" ]; then
            BEST_COUNT=\$COUNT
            BEST_PRIMER=\$pf
        fi
    done

    echo "Best: \$BEST_PRIMER (\$BEST_COUNT pairs)"
    """
}

process REMOVE_PRIMERS {
    tag "${meta.id}"
    cpus 1
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: "*_cutadapt.log", enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/trimmed" : null

    input:
    tuple val(meta), path(r1), path(r2), val(primer_file)

    output:
    tuple val(meta), path("${meta.id}_R1.trimmed.fastq.gz"), path("${meta.id}_R2.trimmed.fastq.gz"), emit: reads
    path("${meta.id}_cutadapt.log"), emit: log

    script:
    """
    cutadapt \\
        -g file:${primer_file} \\
        -G file:${primer_file} \\
        --discard-untrimmed \\
        --pair-filter=any \\
        -j ${task.cpus} \\
        -e ${params.primer_error_rate} \\
        -o ${meta.id}_R1.trimmed.fastq.gz \\
        -p ${meta.id}_R2.trimmed.fastq.gz \\
        ${r1} ${r2} \\
        > ${meta.id}_cutadapt.log 2>&1

    echo "[INFO] ${meta.id}: primer removal complete (${primer_file})" >&2
    """
}
