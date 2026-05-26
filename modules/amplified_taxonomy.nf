nextflow.enable.dsl=2

process amplified_taxonomy {
    label 'short'
    publishDir "${params.output_dir}/report", mode: 'copy'

    input:
        path input_fastas, stageAs: 'input_fastas/*'
        path input_taxonomies, stageAs: 'input_taxonomy/*'
        path homolog_taxids, stageAs: 'homolog_taxids/*'

    output:
        path 'amplified_taxonomy_report.tsv'

    script:
        def ncbi_dir_abs = new File(params.ncbi_dir.toString()).absolutePath

        """
        set -euo pipefail

        ampliphy_amplified_taxonomy.py \
          --ncbi-dir '${ncbi_dir_abs}' \
          --fasta-dir input_fastas \
          --input-taxonomy-dir input_taxonomy \
          --homolog-taxids-dir homolog_taxids \
          --output amplified_taxonomy_report.tsv
        """
}