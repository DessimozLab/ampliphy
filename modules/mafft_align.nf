nextflow.enable.dsl=2

process mafft_align {
    label params.minimal ? 'minimal' : 'short'
    publishDir "${params.output_dir}/msa", mode: 'copy'

    input:
        tuple val(id), path(fa)

    output:
        path "${id}.msa.fa"

    script:
        def mafft_preset = params.mafft_preset.toLowerCase()
        def mafft_options = params.mafft_options
        def threads = params.minimal ? 1 : params.threads

        Map mafft_preset_flags = [
            'auto' : '--auto',
            'fast' : '--retree 2 --maxiterate 0',
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

        if [[ ${fa} == *.gz ]]; then
            gzip -dc ${fa} > "input.fa"
        else
            ln -s ${fa} "input.fa"
        fi

        mafft --quiet --thread ${threads} ${mafft_flag} input.fa \
        | awk '/^>/{if(x) print x; print; x=""; next} {x=x\$0} END{print x}' \
        > ${id}.msa.fa
        """
}
