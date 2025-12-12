nextflow.enable.dsl=2

process mafft_augment {
    label 'short'
    publishDir params.output_dir, mode: 'copy'

    input:
        tuple val(id), path(msa), path(homologs)

    output:
        path "${id}.aug.fa"

    script:
        def mafft_preset  = params.mafft_preset.toLowerCase()
        def mafft_options = params.mafft_options
        def threads       = params.threads ?: 1

        Map mafft_preset_flags = [
            'auto'  : '--auto',
            'fast'  : '--retree 2 --maxiterate 0',
            'linsi' : '--localpair --maxiterate 1000',
            'ginsi' : '--genafpair --maxiterate 1000',
            'einsi' : '--globalpair --maxiterate 1000'
        ]

        if( !mafft_preset_flags.containsKey(mafft_preset) ) {
            error "Invalid MAFFT preset: ${mafft_preset}. Valid options are: ${mafft_preset_flags.keySet().join(', ')}"
        }

        def mafft_flag = (mafft_preset_flags[mafft_preset] + ( mafft_options ? " ${mafft_options}" : "" )).trim()

        """
        set -euo pipefail

        mafft --anysymbol --quiet --keeplength --thread ${threads} ${mafft_flag} \\
          --addfragments "${homologs}" "${msa}" \\
          > ${id}.aug.fa
        """
}
