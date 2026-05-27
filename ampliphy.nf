#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { amplified_taxonomy } from './modules/amplified_taxonomy.nf'
include { homolog_taxonomy } from './modules/homolog_taxonomy.nf'
include { mafft_align } from './modules/mafft_align.nf'
include { mmseqs_prepare_db } from './modules/mmseqs_prepare_db.nf'
include { mmseqs_search } from './modules/mmseqs_search.nf'
include { mafft_amplify } from './modules/mafft_amplify.nf'
include { iqtree_inference } from './modules/iqtree_inference.nf'
include { root_and_prune }  from './modules/root_and_prune.nf'
include { input_lineage } from './modules/input_lineage.nf'
include { amplified_tcs } from './modules/amplified_tcs.nf'

workflow {
    main:
        def input_dir = params.input_dir
        def output_dir = params.output_dir

        new File(output_dir).mkdirs()
        
        def patterns = [
            "${input_dir}/*.{fa,fasta,faa,fna,ffn,frn}{.gz,}",
        ]

        channel
            .fromPath(patterns, checkIfExists: true)
            .ifEmpty { error "No input files found matching: ${input_dir}" }
            .map { file ->
                def name = file.getSimpleName()
                def base = name.replaceFirst(/\.(fa|fasta|fna|ffn|faa|frn)(\.gz)?$/, '')
                tuple( base, file )
            }
            .multiMap { tup -> mafft: tup; mmseqs: tup; taxonomy: tup }
            .set { seq_inputs }

        // log.info "AmpliPhy - MAFFT alignment"
        mafft_align( seq_inputs.mafft )
        // log.info "AmpliPhy - MMseqs2 database preparation"
        def db_channel = mmseqs_prepare_db()
        // log.info "AmpliPhy - MMseqs2 search"
        mmseqs_search(
            seq_inputs.mmseqs,
            db_channel.mmseqs_db_path,
            db_channel.taxonomy_status
        )
        mmseqs_search.out.homolog_search_stats
            .map { report_file -> report_file.text.trim() }
            .collectFile(
                name: 'homolog_search_report.tsv',
                seed: "family_id\toriginal_sequences\thomologs_found\thomologs_added\tamplified_sequences",
                newLine: true,
                sort: true,
                storeDir: "${output_dir}/report"
            )
            
        if( !params.no_taxonomy ) {
            mmseqs_search.out.homolog_taxids
                .collect()
                .set { homolog_taxonomy_inputs }

            homolog_taxonomy( homolog_taxonomy_inputs )
        }
        
        if( params.input_taxonomy ) {
            if( params.no_taxonomy ) {
                error "--input_taxonomy cannot be used with --no_taxonomy because added homolog taxids are required for the amplified taxonomy report."
            }

            channel
                .fromPath("${params.input_taxonomy}/*.tax", checkIfExists: true)
                .ifEmpty { error "No <family_id>.tax files found in --input_taxonomy directory: ${params.input_taxonomy}" }
                .collect()
                .set { input_taxonomy_files }

            seq_inputs.taxonomy
                .map { _id, fasta -> fasta }
                .collect()
                .set { input_fasta_files }

            mmseqs_search.out.homolog_taxids
                .collect()
                .set { added_homolog_taxid_files }

            amplified_taxonomy(
                input_fasta_files,
                input_taxonomy_files,
                added_homolog_taxid_files
            )
        }

        mafft_align.out
            .map { msa_file ->
                def name = msa_file.getSimpleName()
                def base = name.replaceFirst(/\.msa\.fa(\.gz)?$/, '')
                tuple(base, msa_file)
            }
            .set { msa_tuples }

        mmseqs_search.out.homolog_seqs
            .map { hom_file ->
                def name = hom_file.getSimpleName()
                def base = name.replaceFirst(/\.homologs\.fa(\.gz)?$/, '')
                tuple(base, hom_file)
            }
            .set { hom_tuples }

        // Join on id -> (id, msa, homologs)
        msa_tuples
            .join(hom_tuples)
            .set { amplify_inputs }

        // log.info "AmpliPhy - MAFFT amplification with homologs"
        mafft_amplify( amplify_inputs )

        mafft_amplify.out
            .map { amp_file ->
                def name = amp_file.getSimpleName()
                def base = name.replaceFirst(/\.amp\.fa(\.gz)?$/, '')
                def nseq = 0
                amp_file.eachLine { line ->
                    if( line.startsWith('>') ) { nseq = nseq + 1 }
                }
                tuple(base, 'amp.unpruned', amp_file, nseq)
            }
            .set { amp_iqtree_inputs }

        def ori_iqtree_inputs = channel.empty()

        if( params.input_taxonomy ) {
            ori_iqtree_inputs = mafft_align.out.map { msa_file ->
                def name = msa_file.getSimpleName()
                def base = name.replaceFirst(/\.msa\.fa(\.gz)?$/, '')
                def nseq = 0
                msa_file.eachLine { line ->
                    if( line.startsWith('>') ) { nseq = nseq + 1 }
                }
                tuple(base, 'ori', msa_file, nseq)
            }
        }

        amp_iqtree_inputs
            .mix(ori_iqtree_inputs)
            .set { iqtree_inputs }

        iqtree_inference( iqtree_inputs )

        iqtree_inference.out.inferred_trees
            .filter { _id, tree_type, _nwk_file -> tree_type == 'amp.unpruned' }
            .map { id, _tree_type, nwk_file -> tuple(id, nwk_file) }
            .join(hom_tuples)
            .set { prune_inputs }

        root_and_prune( prune_inputs, file(params.mad_script) )

        if( params.input_taxonomy ) {
            input_lineage( input_fasta_files, input_taxonomy_files )

            input_lineage.out.lineages
                .flatten()
                .map { lineage_file ->
                    def name = lineage_file.getSimpleName()
                    def base = name.replaceFirst(/\.lineage\.tsv$/, '')
                    tuple(base, lineage_file)
                }
                .set { lineage_tuples }

            iqtree_inference.out.inferred_trees
                .filter { _id, tree_type, _nwk_file -> tree_type == 'ori' }
                .map { id, _tree_type, nwk_file -> tuple(id, nwk_file) }
                .set { ori_tree_tuples }

            ori_tree_tuples
                .join(root_and_prune.out.amplified_trees)
                .join(lineage_tuples)
                .set { tcs_inputs }

            amplified_tcs( tcs_inputs, file(params.tcs_script) )

            amplified_tcs.out.tcs_rows
                .map { row_file -> row_file.text.trim() }
                .collectFile(
                    name: 'amplified_tcs_report.tsv',
                    seed: "family_id\ttcs_ori\ttcs_amp\n",
                    newLine: true,
                    sort: true,
                    storeDir: "${output_dir}/report"
                )
        }
}