nextflow.enable.dsl=2

process homolog_taxonomy {
    label params.minimal ? 'minimal' : 'short'
    publishDir "${params.output_dir}/report", mode: 'copy'

    input:
        path taxid_files

    output:
        path 'homolog_taxonomy_report.tsv'

    script:
        def ncbi_dir = params.ncbi_dir ?: './ncbi_taxonomy'
        def ncbi_dir_abs = new File(ncbi_dir.toString()).absolutePath

        """
        set -euo pipefail

        ampliphy_homolog_taxonomy.py \
          --ncbi-dir '${ncbi_dir_abs}' \
          --output homolog_taxonomy_report.tsv \
          *.homolog_taxids.tsv
        """
}