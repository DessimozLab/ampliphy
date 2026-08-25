nextflow.enable.dsl=2

process iqtree_inference {
    label params.minimal ? 'minimal' : 'tree'
    publishDir "${params.output_dir}/tree", mode: 'copy',
        saveAs: { filename ->
            params.keep_unpruned_tree && filename.endsWith('.amp.unpruned.nwk') ? filename : null
        }

    input:
        tuple val(id), val(tree_type), path(msa), val(nseq)

    output:
        tuple val(id), val(tree_type), path("${id}.${tree_type}.nwk"), emit: inferred_trees

    script:
        def iqtree_options = params.iqtree_options ?: '-m JTT+I+G4 -B 1000 -keep-ident'
        def threads = params.minimal ? 1 : params.threads

        if( nseq < 4 ) {
            iqtree_options = iqtree_options
                .replaceAll(/(^|\s)(-B|--ufboot|-b|--boot)\s+\S+/, ' ')
                .replaceAll(/(^|\s)-bnni(?=\s|$)/, ' ')
                .replaceAll(/\s+/, ' ')
                .trim()

            log.warn "Bootstrap analysis disabled for ${id}: amplified MSA contains only ${nseq} sequences."
        }

        """
        iqtree -s ${msa} -pre ${id}.${tree_type} -T ${threads} ${iqtree_options}
        cp ${id}.${tree_type}.treefile ${id}.${tree_type}.nwk
        """
}
