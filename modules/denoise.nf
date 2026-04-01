// Amplicon denoising processes.
//
// Three-step workflow:
//   1. FILTER_TRIM    — per-sample quality filtering (uses --lang)
//   2. LEARN_ERRORS   — per-plate error model learning (uses --denoise_engine)
//   3. DENOISE        — per-plate denoising + pair merging (uses --denoise_engine)
//
// --denoise_engine controls which implementation to use:
//   'python' (default) = papa2 (byte-identical to R dada2, no R dependency)
//   'R' = R dada2 package (reference implementation)
//
// --lang controls filter_trim and downstream scripts independently.

// Resolve effective denoise engine: explicit param > lang fallback
def denoiseEngine() { return params.denoise_engine ?: params.lang }

// Output extensions: denoise steps always output .rds when engine=R
def errExt() { return denoiseEngine() == 'python' ? 'pkl' : 'rds' }
def seqExt() { return params.lang == 'python' ? 'pkl' : 'rds' }

process AUTO_TRIM {
    label 'process_low'
    conda params.lang == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml"
    publishDir "${params.outdir}/quality_check", mode: 'copy'

    input:
    path(input_dir)

    output:
    path("auto_trim.tsv"), emit: params_tsv
    env(TRUNC_FWD), emit: trunc_fwd
    env(TRUNC_REV), emit: trunc_rev

    script:
    """
    microscape auto-trim "${input_dir}" \
        --min-quality ${params.auto_trim_min_quality} \
        --output auto_trim.tsv \
        --verbose

    export TRUNC_FWD=\$(grep trunc_len_fwd auto_trim.tsv | cut -f2)
    export TRUNC_REV=\$(grep trunc_len_rev auto_trim.tsv | cut -f2)
    echo "[INFO] Auto-trim: fwd=\$TRUNC_FWD rev=\$TRUNC_REV"
    """
}

process FILTER_TRIM {
    tag "${meta.id}"
    label 'process_low'
    conda params.lang == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml"
    publishDir "${params.outdir}/filtered", mode: 'copy', pattern: "*_filt_stats.tsv", enabled: !params.store_dir

    input:
    tuple val(meta), path(r1), path(r2)
    val(trunc_fwd)
    val(trunc_rev)

    output:
    tuple val(meta), path("${meta.id}_R1.filt.fastq.gz"), path("${meta.id}_R2.filt.fastq.gz"), emit: reads
    path("${meta.id}_filt_stats.tsv"), emit: stats

    script:
    def eff_trunc_fwd = params.truncLen_fwd > 0 ? params.truncLen_fwd : trunc_fwd
    def eff_trunc_rev = params.truncLen_rev > 0 ? params.truncLen_rev : trunc_rev
    if (params.lang == 'python')
    """
    papa2 filter-trim \
        "${r1}" "${meta.id}_R1.filt.fastq.gz" \
        "${r2}" "${meta.id}_R2.filt.fastq.gz" \
        --max-ee ${params.maxEE} \
        --trunc-q ${params.truncQ} \
        --max-n ${params.maxN} \
        --trunc-len-fwd ${eff_trunc_fwd} \
        --trunc-len-rev ${eff_trunc_rev} \
        --stats "${meta.id}_filt_stats.tsv" \
        --sample-id "${meta.id}"
    """
    else
    """
    dada2_filter_trim.R \
        "${meta.id}" "${r1}" "${r2}" \
        ${params.maxEE} ${params.truncQ} ${params.maxN} \
        ${params.truncLen_fwd} ${params.truncLen_rev} \
        ${task.cpus}
    """
}

// Learn error rates — uses denoise_engine (R or python/papa2)
process LEARN_ERRORS {
    tag "${meta.id}"
    label 'process_high'
    conda ((params.denoise_engine ?: params.lang) == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml")
    publishDir "${params.outdir}/error_models", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/error_models" : null

    input:
    tuple val(meta), path(r1_files), path(r2_files)

    output:
    tuple val(meta.plate), val(meta), path("${meta.id}_errF.${errExt()}"), path("${meta.id}_errR.${errExt()}"), emit: error_models
    path("${meta.id}_error_rates.pdf"), emit: error_plots

    script:
    if (denoiseEngine() == 'python')
    """
    learn_errors.py "${meta.id}" ${task.cpus}
    """
    else
    """
    dada2_learn_errors.R "${meta.id}" ${task.cpus}
    """
}

// Per-plate denoising — uses denoise_engine (R or python/papa2)
// Output is always .pkl when lang=python (converted from .rds by wrapper if needed)
process DENOISE {
    tag "${meta.id}"
    label 'process_high'
    conda ((params.denoise_engine ?: params.lang) == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml")
    publishDir "${params.outdir}/seqtabs", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/seqtabs" : null

    input:
    tuple val(meta), path(r1_files), path(r2_files), path(errF), path(errR)

    output:
    tuple val(meta), path("${meta.id}.seqtab.${seqExt()}"), emit: seqtab
    path("${meta.id}.seqtab.tsv"), emit: seqtab_tsv

    script:
    if (denoiseEngine() == 'python')
    """
    denoise.py \
        "${meta.id}" "${errF}" "${errR}" \
        ${params.min_overlap} ${task.cpus}
    """
    else if (params.lang == 'python')
    // R dada2 engine but Python downstream: run R, then convert .rds to .pkl
    """
    dada2_denoise.R \
        "${meta.id}" "${errF}" "${errR}" \
        ${params.min_overlap} ${task.cpus}

    rds_to_pkl.py "${meta.id}.seqtab.rds" "${meta.id}.seqtab.pkl"
    """
    else
    """
    dada2_denoise.R \
        "${meta.id}" "${errF}" "${errR}" \
        ${params.min_overlap} ${task.cpus}
    """
}
