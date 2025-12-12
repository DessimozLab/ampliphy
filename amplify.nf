#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { mafft_align } from './modules/mafft_align.nf'
include { mmseqs_prepare_db } from './modules/mmseqs_prepare_db.nf'
include { mmseqs_search } from './modules/mmseqs_search.nf'
include { mafft_augment } from './modules/mafft_augment.nf'

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
            .set { augment_inputs }

        // log.info "AmpliPhy - MAFFT augmentation with homologs"
        mafft_augment( augment_inputs )
}