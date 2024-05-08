#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../modules/nf-core/diamond/makedb/main' 
include { BLAST_UPDATEBLASTDB } from '../modules/nf-core/blast/updateblastdb'
include { BLAST_BLASTDBCMD } from '../modules/nf-core/blast/blastdbcmd'
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
    DOWNLOAD_OBOFILE()
    BLAST_UPDATEBLASTDB()
    DIAMOND_MAKEDB()
    // TODO: Download InterProScan data
    DOWNLOAD_KO()
}

// TODO: Move download stuff into modules

process oboFile {

  publishDir params.oboFolder, mode: 'copy'
  label 'download'

  output:
  file "gene_ontology.obo" into oboFile

  """
  curl --retry 3 -o gene_ontology.obo ${params.goOboURL};
  """

}


process downloadInterPro {

  publishDir params.dbipscanPath, mode: 'copy'
  label 'ipscan'

  output:
  file "*" into interpro_data

  """
  curl --retry 3 -o iprscan.tar.gz ${params.iprscanURL};
  tar zxf iprscan.tar.gz
  rm iprscan.tar.gz
  cd interproscan-${params.iprscanVersion}
  python3 initial_setup.py
  cd ..
  mv interproscan-${params.iprscanVersion}/data .
  rm -rf interproscan-${params.iprscanVersion}
  mv data/* .
  """

}

process downloadKO {

  publishDir params.dbKOPath, mode: 'copy'

  label 'download'

  output:
  file "ko_list" into ko_list
  file "profiles" into ko_profiles
  file "ko_store" into ko_store

  """
  curl --retry 3 -o ko_list.gz ${params.koURLlist};
  gunzip ko_list.gz;
  curl --retry 3 -o profiles.tar.gz ${params.koURLprofiles};
  tar zxf profiles.tar.gz; rm profiles.tar.gz;
  mkdir ko_store
  bulkDownloadKEGG.pl ko_list ko_store
  """

}

