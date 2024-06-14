#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../modules/nf-core/diamond/makedb/main'
include { BLAST_UPDATEBLASTDB } from '../modules/nf-core/blast/updateblastdb'
include { BLAST_BLASTDBCMD } from '../modules/nf-core/blast/blastdbcmd'
include { INTERPROSCAN_DOWNLOAD } from '../modules/local/fa-nf/interproscan/download'
include { KOFAMSCAN_DOWNLOAD } from '../modules/local/fa-nf/kofamscan/download'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_fanf_pipeline'

if ( params.dbPath == null || params.dbPath == "" ) {
  log.info "No target directory specified"
  exit 1
}

if ( params.dbList == null || params.dbList == "" ) {
  log.info "No BLAST DBs provided"
  exit 1
}

workflow DOWNLOAD {

    take:
    dbnames // channel: [ dbnames ]

    main:
    ch_versions = Channel.empty()

    // DOWNLOAD_OBOFILE()

    // Let's map to fit into module requirement
    // TODO: Need to check why trimming is required here
    dbmap = dbnames.map { it -> [ [ id: it.trim() ], it.trim() ] }
    BLAST_UPDATEBLASTDB(dbmap) // We should convert this channel of ids into meta
    // BLAST_UPDATEBLASTDB.out.db.view()

    blastmap =  BLAST_UPDATEBLASTDB.out.db.map{ it -> [ it[0], 'all', [] ] }
    // blastmap.view()
    BLAST_BLASTDBCMD( blastmap, BLAST_UPDATEBLASTDB.out.db )
    DIAMOND_MAKEDB(BLAST_BLASTDBCMD.out.fasta, [], [], [])
    // TODO: Download InterProScan data
    // DOWNLOAD_INTERPROSCAN()
    // KOFAMSCAN_DOWNLOAD()

    ch_versions = ch_versions.mix(BLAST_UPDATEBLASTDB.out.versions.first())
    ch_versions = ch_versions.mix(BLAST_BLASTDBCMD.out.versions.first())

    emit:
    fasta = BLAST_BLASTDBCMD.out.fasta
    fasta.view { println "Path: ${it}" }
    versions = ch_versions
    versions.view { println "Versions: ${it}" }
}

// TODO: Move download stuff into modules

// process oboFile {
//
//   publishDir params.oboFolder, mode: 'copy'
//   label 'download'
//
//   output:
//   file "gene_ontology.obo" into oboFile
//
//   """
//   curl --retry 3 -o gene_ontology.obo ${params.goOboURL};
//   """
//
// }
//
