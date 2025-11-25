#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { mafft_align } from './modules/mafft_align.nf'

workflow {
    main:
        def input_dir = params.input_dir
        def output_dir = params.output_dir

        new File(output_dir).mkdirs()
        log.info "AmpliPhy - MAFFT alignment"

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
            .set { seq_files }

        mafft_align( seq_files )
}