nextflow.enable.dsl=2

process input_lineage {
    label 'short'

    input:
        path fasta_files, stageAs: 'fastas/*'
        path taxonomy_files, stageAs: 'taxonomy/*'

    output:
        path 'lineages/*.lineage.tsv', emit: lineages

    script:
        def ncbi_dir_abs = new File(params.ncbi_dir.toString()).absolutePath

        """
        ampliphy_input_lineage.py --ncbi-dir '${ncbi_dir_abs}' --fasta-dir fastas --taxonomy-dir taxonomy --output-dir lineages
        """
}