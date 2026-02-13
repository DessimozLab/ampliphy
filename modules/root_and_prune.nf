nextflow.enable.dsl=2

process root_and_prune {
    label params.minimal ? 'minimal' : 'short'
    publishDir "${params.output_dir}/tree", mode: 'copy'

    input:
        tuple val(id), path(amp_nwk), path(homologs_fa)
        path mad_py

    output:
        path "${id}.nwk"

    script:
        def no_rooting = params.no_rooting ? 'true' : 'false'

        """
        set -euo pipefail

        ROOTED="${id}.amp.nwk.rooted"

        if [[ "${no_rooting}" == "true" ]]; then
            cp "${amp_nwk}" "\${ROOTED}"
        else
            "./${mad_py}" "${amp_nwk}"
        fi

        TIPFILE="${id}.tips.txt"
        awk '/^>/{print substr(\$1,2)}' "${homologs_fa}" > "\${TIPFILE}"

        gotree prune -i "\${ROOTED}" -f "\${TIPFILE}" -o "${id}.nwk"
        """
}
