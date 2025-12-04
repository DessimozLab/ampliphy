nextflow.enable.dsl=2

process mmseqs_search {
    label 'mmseqs'
    publishDir params.output_dir, mode: 'copy'

    input:
        tuple val(id), path(fa)

    output:
        path "${id}.homologs.fa", emit: homolog_seqs

    script:
        def database = (params.database ?: 'UniRef50').toString()
        def custom_database = params.custom_database ?: ''
        def tmp_root = params.tmp_dir ?: './tmp'
        def db_root = params.database_dir ?: './mmseqs_db'
        def mmseqs_options = params.mmseqs_options ?: ''
        def max_depth_param = (params.max_depth ?: 5).toString()
        def max_seqs_param = (params.max_seqs ?: 0).toString()
        def threads = params.threads ?: 4

        def allowed_dbs = ['UniRef100', 'UniRef90', 'UniRef50', 'UniProtKB', 'UniProtKB/Swiss-Prot', 'NR']

        if( !custom_database && !allowed_dbs.contains(database) ) {
            error "Invalid MMseqs2 database: ${database}. Valid options are: ${allowed_dbs.join(', ')} or provide a custom database path via --custom_database"
        }

        // Prepare paths
        def custom_database_absolute = custom_database ? new File(custom_database.toString()).absolutePath : ''
        def db_root_absolute = new File(db_root.toString()).absolutePath
        def tmp_root_absolute = new File(tmp_root.toString()).absolutePath
        def db_path = custom_database_absolute ?: "${db_root_absolute}/${database.replaceAll(/[\\\\/:\\s]+/, '_')}"

        """
        set -euo pipefail

        ampliphy_mmseqs_search.sh \\
            --id ${id} \\
            --input ${fa} \\
            --database '${database}' \\
            --db-path '${db_path}' \\
            --custom-db '${custom_database_absolute}' \\
            --tmp-root '${tmp_root_absolute}' \\
            --mmseqs-options '${mmseqs_options}' \\
            --max-depth ${max_depth_param} \\
            --max-seqs ${max_seqs_param} \\
            --threads ${threads}
        """
}