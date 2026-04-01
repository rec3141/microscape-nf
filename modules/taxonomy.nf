// Taxonomy assignment — naive Bayesian kmer classifier.
//
// Each reference database runs as an independent task, enabling parallel
// classification across databases. Supports --lang R (dada2) or --lang python
// (custom Wang et al. 2007 implementation with multiprocessing).

def taxExt() { return params.lang == 'python' ? 'pkl' : 'rds' }

process ASSIGN_TAXONOMY {
    tag "${db_name}"
    label 'process_high'
    conda params.lang == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml"
    publishDir "${params.outdir}/taxonomy", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/taxonomy" : null

    input:
    tuple val(db_name), path(ref_db), val(tax_levels)
    path(seqtab)

    output:
    tuple val(db_name), path("${db_name}_taxonomy.${taxExt()}"), path("${db_name}_bootstrap.${taxExt()}"), emit: taxonomy
    path("${db_name}_taxonomy.tsv"), emit: taxonomy_tsv

    script:
    def levels_arg = tax_levels ? "\"${tax_levels}\"" : "null"
    def levels_py  = tax_levels ? "'" + tax_levels.split(",").collect{ it.trim() }.join("','") + "'" : "None"
    if (params.lang == 'python')
    """
    python3 -c "
import os, pickle, pandas as pd
os.environ['OMP_NUM_THREADS'] = '${task.cpus}'
from papa2.taxonomy import assign_taxonomy

data = pickle.load(open('${seqtab}', 'rb'))
seqs = sorted(data['sequence'].unique()) if isinstance(data, pd.DataFrame) else list(data)
print(f'[INFO] {len(seqs)} query sequences')

levels = [${levels_py}] if ${levels_py} is not None else None
kw = dict(seqs=seqs, ref_fasta='${ref_db}', min_boot=50, output_bootstraps=True, verbose=True)
if levels: kw['tax_levels'] = levels
result = assign_taxonomy(**kw)

tax_df, boot_df = result['tax'], result['boot']
pickle.dump(result, open('${db_name}_taxonomy.pkl', 'wb'))
pickle.dump(boot_df, open('${db_name}_bootstrap.pkl', 'wb'))
tax_out = tax_df.copy(); tax_out.insert(0, 'sequence', tax_out.index)
tax_out.to_csv('${db_name}_taxonomy.tsv', sep='\t', index=False)
print(f'[INFO] Taxonomy complete for ${db_name}')
"
    """
    else
    """
    assign_taxonomy.R \
        "${seqtab}" "${ref_db}" "${db_name}" \
        ${task.cpus} ${levels_arg}
    """
}
