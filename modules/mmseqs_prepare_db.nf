nextflow.enable.dsl=2

process mmseqs_prepare_db {
    label params.minimal ? 'minimal' : 'mmseqs'

    output:
        path 'mmseqs_db_target.txt', emit: mmseqs_db_path

    script:
        def custom_database = params.custom_database ?: ''
        def tmp_root        = params.tmp_dir ?: './tmp'
        def db_root         = params.database_dir ?: './mmseqs_db'

        // Raw user input (for messages)
        def db_raw   = (params.database ?: 'UniRef50').toString().trim()
        def db_lower = db_raw.toLowerCase()

        // Case-insensitive alias map -> canonical MMseqs DB name
        def db_aliases = [
            'uniref100'              : 'UniRef100',
            'uniref90'               : 'UniRef90',
            'uniref50'               : 'UniRef50',
            'uniprotkb'              : 'UniProtKB',
            'uniprotkb/swiss-prot'   : 'UniProtKB/Swiss-Prot',
            'uniprotkb_swiss-prot'   : 'UniProtKB/Swiss-Prot',
            'uniprotkb swiss-prot'   : 'UniProtKB/Swiss-Prot',
            'swiss-prot'             : 'UniProtKB/Swiss-Prot',
            'swissprot'              : 'UniProtKB/Swiss-Prot',
            'nr'                     : 'NR',
        ]
        
        // Canonical MMseqs database name (what we pass to `mmseqs databases`)
        def mmseqs_db_name = db_aliases[db_lower]

        if( !custom_database && !mmseqs_db_name ) {
        def valid = [
            'UniRef100',
            'UniRef90',
            'UniRef50',
            'UniProtKB',
            'UniProtKB/Swiss-Prot', 'Swiss-Prot', 'SwissProt',
            'NR'
        ]
        error "Invalid MMseqs2 database: '${db_raw}'. " +
            "Valid options are (case-insensitive): " +
            valid.join(', ') +
            ", or provide a custom database path via --custom_database. "
        }

        def custom_database_absolute = custom_database ? new File(custom_database.toString()).absolutePath : ''
        def db_root_abs  = new File(db_root.toString()).absolutePath
        def tmp_root_abs = new File(tmp_root.toString()).absolutePath

        def db_name_for_mmseqs = mmseqs_db_name ?: db_raw
        def db_dir_name = mmseqs_db_name
            ? mmseqs_db_name.replaceAll(/[^A-Za-z0-9_.-]+/, '_')
            : 'custom_db'

        // Final DB prefix on disk
        def db_target = custom_database_absolute ?: "${db_root_abs}/${db_dir_name}"

        """
        set -euo pipefail

        ampliphy_mmseqs_prepare_db.sh \
        --database '${db_name_for_mmseqs}' \
        --db-target '${db_target}' \
        --custom-db '${custom_database_absolute}' \
        --tmp-root '${tmp_root_abs}'

        printf '%s\\n' '${db_target}' > mmseqs_db_target.txt
        """
}
