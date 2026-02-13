nextflow.enable.dsl=2

process mmseqs_search {
    label params.minimal ? 'minimal' : 'mmseqs'
    publishDir "${params.output_dir}/homologs", mode: 'copy'

    input:
        tuple val(id), path(fa)
        path db_file

    output:
        path "${id}.homologs.fa", emit: homolog_seqs

    script:
        def tmp_root        = params.tmp_dir ?: './tmp'
        def mmseqs_options  = params.mmseqs_options ?: ''
        def max_depth_param = (params.max_depth ?: 5).toString()
        def max_seqs_param  = (params.max_seqs ?: 0).toString()
        def threads         = params.minimal ? 1 : params.threads

        def tmp_root_abs = new File(tmp_root.toString()).absolutePath

        """
        set -euo pipefail

        DB_TARGET=\$(cat ${db_file})

        ampliphy_mmseqs_search.sh \
          --id ${id} \
          --input ${fa} \
          --db-target "\${DB_TARGET}" \
          --tmp-root '${tmp_root_abs}' \
          --mmseqs-options '${mmseqs_options}' \
          --max-depth '${max_depth_param}' \
          --max-seqs '${max_seqs_param}' \
          --threads ${threads}
        """
}
