nextflow.enable.dsl=2

process iqtree_inference {
    label params.minimal ? 'minimal' : 'short'
    publishDir "${params.output_dir}/tree", mode: 'copy'

    input:
        tuple val(id), path(amp_fa)

    output:
        path "${id}.amp.nwk"

    script:
        def iqtree_options = params.iqtree_options ?: '-m JTT+I+G4 -B 1000'
        def threads = params.minimal ? 1 : params.threads

        """
        set -euo pipefail
        iqtree -s "${amp_fa}" -pre "${id}.amp" -T ${threads} ${iqtree_options}
        cp "${id}.amp.treefile" "${id}.amp.nwk"
        """
}
