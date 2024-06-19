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
     .splitFasta( by: 10, file: true )
     .map  { it -> [ [ id: it.toString().trim().split("/")[-1] ], it ] }

    ch_fasta.view()

    if ( ! params.outdir ) {
        log.error "No DIR"
        exit 1
    }

    ch_db = Channel
    .fromPath( "${params.outdir}/diamond/${params.dbname}.dmnd", checkIfExists:true )
     .map  { it -> [ [ id: it.toString().trim().split("/")[-1] ], it ] }

    ch_diamond = DIAMOND_BLASTP(ch_fasta, ch_db, "blast", '')
    ch_ko_profiles = Channel.fromPath("$params.outdir/kofamscan/profiles", checkIfExists: true)
    ch_ko_list = Channel.fromPath("$params.outdir/kofamscan/ko_list", checkIfExists: true)

    ch_kofamscan = KOFAMSCAN(ch_fasta, ch_ko_profiles, ch_ko_list)

    ch_versions = Channel.empty()

    emit:
    versions = ch_versions
    ch_fasta
}

