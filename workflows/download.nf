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

if ( params.blastDBList == null || params.blastDBList == "" ) {
  log.info "No BLAST DBs provided"
  exit 1
}

blastDBChannel = Channel.fromList( params.blastDBList?.tokenize(',') )

workflow DOWNLOAD {
    // DOWNLOAD_OBOFILE()
    meta_map = [
        [id: 'mito'], // meta map
        'mito'
    ]
    BLAST_UPDATEBLASTDB(meta_map)
    // DIAMOND_MAKEDB()
    // TODO: Download InterProScan data
    // DOWNLOAD_INTERPROSCAN()
    // KOFAMSCAN_DOWNLOAD()
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
