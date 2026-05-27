nextflow.enable.dsl=2

process amplified_tcs {
    label 'short'

    input:
        tuple val(id), path(ori_tree), path(amp_tree), path(lineage_file)
        path tcs_script

    output:
        path "${id}.tcs.tsv", emit: tcs_rows

    script:
        """
        TCS_ORI=\$(python3 ${tcs_script} ${ori_tree} ${lineage_file} | cut -f2)
        TCS_AMP=\$(python3 ${tcs_script} ${amp_tree} ${lineage_file} | cut -f2)
        printf '%s\\t%s\\t%s\\n' '${id}' "\${TCS_ORI}" "\${TCS_AMP}" > ${id}.tcs.tsv
        """
}