#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { DIAMOND_MAKEDB } from '../modules/nf-core/diamond/makedb/main' 
include { BLAST_UPDATEBLASTDB } from '../modules/nf-core/blast/updateblastdb'
include { BLAST_BLASTDBCMD } from '../modules/nf-core/blast/blastdbcmd'

// default parameters
params.help = false

// Main result and log dirs
params.dbPath = "/nfs/db"

// Version
params.iprscanVersion = "5.52-86.0"
params.koVersion = "2021-05-02"

// URLs
params.iprscanURL = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${params.iprscanVersion}/interproscan-${params.iprscanVersion}-64-bit.tar.gz "
params.koURLlist = "ftp://ftp.genome.jp/pub/db/kofam/archives/${params.koVersion}/ko_list.gz"
params.koURLprofiles = "ftp://ftp.genome.jp/pub/db/kofam/archives/${params.koVersion}/profiles.tar.gz"
params.goOboURL = "http://www.geneontology.org/ontology/gene_ontology.obo"

// NCBI DB list
params.blastDBList = "swissprot,pdbaa"
params.blastTimeout = 600

// Specific DB Paths
Date date = new Date()
String datePart = date.format("yyyyMM")
params.blastDbFolder = "${params.dbPath}/ncbi/${datePart}/blastdb/db"
params.dbipscanPath = "${params.dbPath}/iprscan/${params.iprscanVersion}"
params.dbKOPath = "${params.dbPath}/kegg/${params.koVersion}"
params.oboFolder = "${params.dbPath}/geneontology/${datePart}"

// Mail for sending reports
params.email = ""

//print usage
if ( params.help ) {
  log.info ''
  log.info 'Functional Annotation - Download datasets pipeline'
  log.info '----------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info "  ./nextflow run download.nf --config params.config [options]"
  log.info ''
  log.info 'Options:'
  log.info '-resume		resume pipeline from the previous step, i.e. in case of error'
  log.info '-help		this message'
  exit 1
}

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
  // Download GO OBO file
  // call oboFile
  // // Download NCBI DBs
  // call downloadNCBI
  // // Format NCBI DBs
    DIAMOND_MAKEDB()
  // // Download InterProScan data
  // call downloadInterPro
  // // Download KEGG Orthology data
  // call downloadKO
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


// process downloadNCBI {
//
//   publishDir params.blastDbFolder, mode: 'copy'
//
//   label 'blast'
//
//   input:
//   val db from blastDBChannel
//
//   output:
//   set val(db), file ("*") into blastdb
//
//   """
//   update_blastdb.pl ${db} --timeout ${params.blastTimeout} --decompress
//   blastdbcmd -dbtype prot -db ${db} -entry all -out ${db}.fa
//   """
//
// }
//
// process formatDIAMOND {
//
//   publishDir params.blastDbFolder, mode: 'copy'
//
//   label 'diamond'
//
//   input:
//   set db, file(fasta) from blastdb
//
//   output:
//   file "*" into formatted_blastdb
//
//   """
//   diamond makedb --in ${db}.fa --db ${db}
//   """
//
// }

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


// On finising
workflow.onComplete {

    println ( workflow.success ? "\nDone! Check downloaded datasets in --> $params.dbPath\n" : "Oops .. something went wrong" )

    def msg = """\
    Pipeline execution summary
    ---------------------------
    FA-nf Version    : ${workflow.manifest.version}
    Nextflow Version : ${nextflow.version}
    Command Line     : ${workflow.commandLine}
    Resumed          : ${workflow.resume}
    Completed At     : ${workflow.complete}
    Duration         : ${workflow.duration}
    Success          : ${workflow.success}
    Exit Code        : ${workflow.exitStatus}
    Error Report     : ${workflow.errorReport ?: '-'}
    Launch Dir       : ${workflow.launchDir}
    """
    .stripIndent()

    println( msg )

    if ( mysql ) {
        def procfile = new File( params.mysqllog+"/PROCESS" )
        procfile.delete()
    }

    if (params.email == "yourmail@yourdomain" || params.email == "") {
        log.info 'Skipping email\n'
    } else {
        log.info "Sending email to ${params.email}\n"
        sendMail(to: params.email, subject: "[FA-nf] Execution finished", body: msg)
    }

}

workflow.onError {

 println( "Something went wrong" )

}



