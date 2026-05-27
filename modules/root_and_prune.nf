nextflow.enable.dsl=2

process root_and_prune {
    label 'short'
    publishDir "${params.output_dir}/tree", mode: 'copy'

    input:
        tuple val(id), path(amp_unpruned_nwk), path(homologs_fa)
        path mad_script

    output:
        tuple val(id), path("${id}.amp.nwk"), emit: amplified_trees

    script:
        def no_rooting = params.no_rooting ? 'true' : 'false'

        """
        set -euo pipefail

        ROOTED="${id}.amp.unpruned.nwk.rooted"

        if [[ "${no_rooting}" == "true" ]]; then
            cp "${amp_unpruned_nwk}" "\${ROOTED}"
        else
            chmod +x "${mad_script}" || true
            "./${mad_script}" "${amp_unpruned_nwk}" -n
        fi

        if [[ ! -s "${homologs_fa}" ]]; then
            cp "\${ROOTED}" "${id}.amp.nwk"
        else
            awk '/^>/{print substr(\$1,2)}' "${homologs_fa}" > "${id}.tips.txt"
            gotree prune -i "\${ROOTED}" -f "${id}.tips.txt" -o "${id}.amp.nwk"
        fi
        """
}