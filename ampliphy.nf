#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { mafft_align } from './modules/mafft_align.nf'
include { mmseqs_prepare_db } from './modules/mmseqs_prepare_db.nf'
include { mmseqs_search } from './modules/mmseqs_search.nf'
include { mafft_amplify } from './modules/mafft_amplify.nf'
include { iqtree_inference } from './modules/iqtree_inference.nf'
include { root_and_prune }  from './modules/root_and_prune.nf'

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
            .multiMap { tup -> mafft: tup; mmseqs: tup }
            .set { seq_inputs }

        // log.info "AmpliPhy - MAFFT alignment"
        mafft_align( seq_inputs.mafft )
        // log.info "AmpliPhy - MMseqs2 database preparation"
        def db_channel = mmseqs_prepare_db()
        // log.info "AmpliPhy - MMseqs2 search"
        mmseqs_search( seq_inputs.mmseqs, db_channel.mmseqs_db_path )
        mmseqs_search.out.homolog_search_stats
            .map { report_file -> report_file.text.trim() }
            .collectFile(
                name: 'homolog_search_report.tsv',
                seed: "family_id\toriginal_sequences\thomologs_found\thomologs_added\tamplified_sequences",
                newLine: true,
                sort: true,
                storeDir: "${output_dir}/report"
            )

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
                    if( line.startsWith('>') ) nseq = nseq + 1
                }

                tuple(base, amp_file, nseq)
            }
            .set { amp_tuples }
            
        // log.info "AmpliPhy - IQ-TREE2 phylogenetic inference"    
        iqtree_inference( amp_tuples )

        iqtree_inference.out
            .map { nwk_file ->
            def name = nwk_file.getSimpleName()
                def base = name.replaceFirst(/\.amp\.nwk(\.gz)?$/, '')
                tuple(base, nwk_file)
            }
            .set { tree_tuples }

        // log.info "AmpliPhy - Rooting and pruning trees"
        tree_tuples
            .join(hom_tuples)
            .set { prune_inputs }
        
        root_and_prune( prune_inputs, file(params.mad_script) )
}