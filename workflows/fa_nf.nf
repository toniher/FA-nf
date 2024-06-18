#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// TODO: Validate inputs here

//  Include modules below
include { AGAT_CONVERTSPGXF2GXF } from '../modules/nf-core/agat/convertspgxf2gxf/main'
include { DIAMOND_BLASTP } from '../modules/nf-core/diamond/blastp/main'
include { KOFAMSCAN } from '../modules/nf-core/kofamscan/main'
include { INTERPROSCAN } from '../modules/nf-core/interproscan/main'

workflow FA_NF {

    take:
    fasta

    main:
    ch_fasta = Channel
     .fromPath(fasta, checkIfExists:true)
     .splitFasta( record: [id: true, seqString: true], by: 10 )

    ch_fasta.view()

    if ( ! params.outdir ) {
        log.error "No DIR"
        exit 1
    }

    ch_db = Channel.fromPath( "${params.outdir}/diamond/${params.dbname}.dmnd", checkIfExists:true )

    ch_diamond = DIAMOND_BLASTP(ch_fasta, ch_db, "blast", '')

    ch_versions = Channel.empty()

    emit:
    versions = ch_versions
    ch_fasta
}

