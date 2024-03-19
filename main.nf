#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017-2024, Centre for Genomic Regulation (CRG)
 *
 * Copyright (c) 2017, Anna Vlasova
 *
 * Copyright (c) 2017, Emilio Palumbo
 *
 * Copyright (c) 2018-2024, Toni Hermoso Pulido
 *
 * Functional Annotation Pipeline for protein annotation from non-model organisms
 * from Genome Annotation Team in Catalonia (GATC) implemented in Nextflow
 *
 */

// DSL2 pipeline
nextflow.enable.dsl = 2

include { FA-NF } from './workflows/fa-nf'
include { DOWNLOAD } from './workflows/download'


// default parameters
params.help = false
params.debug = false
params.dbEngine = "SQLite" // MySQL otherwise
params.mysqllog = "${baseDir}/tmp"

// Main input files
params.proteinFile = null;
params.gffFile = null;
params.speciesName = "organism"
params.dbname = "organismDB"

// Main result and log dirs
params.resultPath = "${baseDir}/results/"
params.stdoutLog = "${baseDir}/logs/out.log"
params.stderrLog = "${baseDir}/logs/err.log"

// Sizes for different programs
params.chunkIPSSize = null
params.chunkBlastSize = null
params.chunkKoalaSize = null
params.chunkWebSize = null
params.debugSize = 2

// Blast related params
params.blastFile = null
params.evalue = 0.00001
params.diamond = null
params.blastDbPath = null

// GO retrieval params
params.gogourl = ""
params.gogohits = 30
params.blastAnnotMode = "common" // common, most, all available so far

// KEGG
params.kolist = null
params.koprofiles = null
params.koentries = null
params.kegg_release = null
params.kegg_species = "hsa, dme, cel, ath"
params.keggFile = null

// Params for InterProScan
//  Temporary location for InterproScan intermediary files. This can be huge
params.ipscantmp = "${baseDir}/tmp/"
//  Location of InterproScan properties. Do not modify unless it matches your container
params.ipscanproperties = "/usr/local/interproscan/interproscan.properties"

params.ipscandata = ""
ipscandata = "/usr/local/interproscan/data"
if ( params.ipscandata != "" ) {
  ipscandata = params.ipscandata
} else {
  if ( params.iprscanVersion ) {
    if ( params.dbPath ) {
      if ( new File( params.dbPath + "/iprscan/" + params.iprscanVersion ).exists() ) {
        ipscandata = "${params.dbPath}/iprscan/${params.iprscanVersion}"
      }
    }
  }
}

// Params for dealing with GFF
params.gffclean = true
params.gffstats = true
// Remove version from protein entries (e.g. X5543AP.2)
params.rmversion = false

// File with GO information, otherwise is downloaded
params.oboFile = null

// Skip params
params.skip_cdSearch = false
params.skip_sigtarp = false

// Mail for sending reports
params.email = ""

// Download paths
params.dbPath = null
params.blastDbFolder = null
params.dbipscanPath = null
params.dbKOPath = null

// Handling defaults from download part
Date date = new Date()
String datePart = date.format("yyyyMM")

// Handle BlastDB paths
blastDbPath = null;

if ( ! params.blastDbPath ) {

  if ( ! params.blastDbFolder ) {

    if ( ! params.dbPath ) {

      log.info "At least dbPath needs to be defined!"
      exit 1

    } else {
      blastDbFolder = "${params.dbPath}/ncbi/${datePart}/blastdb/db"
      // Let's assume Swissprot
      blastDbPath = "${blastDbFolder}/swissprot"
    }

  } else {
    // Let's assume Swissprot
    blastDbPath = "${params.blastDbFolder}/swissprot"
  }

} else {
  blastDbPath = params.blastDbPath
}

// KEGG paths
kolist = null
koprofiles = null
koentries = null

if ( params.kolist ) {
  kolist = params.kolist
}

if ( params.koprofiles ) {
  koprofiles = params.koprofiles
}

if ( params.koentries ) {
  koentries = params.koentries
}

if ( ! params.kolist && ! params.keggFile ) {
  if ( params.koVersion && params.dbPath ) {
    kolist = "${params.dbPath}/kegg/${params.koVersion}/ko_list"
  }
}

if ( ! params.koprofiles && ! params.keggFile ) {
  if ( params.koVersion && params.dbPath ) {
    koprofiles = "${params.dbPath}/kegg/${params.koVersion}/profiles"
  }
}

if ( ! params.koentries ) {
  if ( params.koVersion && params.dbPath ) {
    koentries = "${params.dbPath}/kegg/${params.koVersion}/ko_store"
  }
}

//print usage
if ( params.help ) {
  log.info ''
  log.info 'Functional annotation pipeline'
  log.info '----------------------------------------------------'
  log.info 'Run functional annotation for a given species.'
  log.info ''
  log.info 'Usage: '
  log.info "  ./nextflow run main.nf --config params.config [options]"
  log.info ''
  log.info 'Options:'
  log.info '-resume		resume pipeline from the previous step, i.e. in case of error'
  log.info '-help		this message'
  exit 1
}

/*
* Parse the input parameters
*/

// species-specific parameters
protein = null
annotation = null
config_file = file(params.config)

dbFile = false
boolean exists = false
boolean mysql = false
gffavail = false
gffclean = false
gffstats = false
// Skip
skip_cdSearch = false
skip_sigtarp = false

if( params.dbEngine.toLowerCase()=="mysql" ) {
 mysql = true
}

if ( params.gffclean ) {
 gffclean = true
}

if ( params.gffstats ) {
 gffstats = true
}

if ( params.skip_cdSearch ) {
 skip_cdSearch = true
}

if ( params.skip_sigtarp ) {
 skip_sigtarp = true
}

// Handling MySQL in a cleaner way
dbhost = null

// Getting contents of file
if ( mysql ) {
 dbhost = "127.0.0.1" // Default value. Localhost

 if ( new File(  params.mysqllog+"/DBHOST" ).exists() ) {
  dbhost = new File(  params.mysqllog+"/DBHOST" ).text.trim()
 }
} else {
 dbFileName = params.resultPath+params.dbname+'.db'
 dbFile = file(dbFileName)
 if ( dbFile.exists() && dbFile.size() > 0 ) {
  exists = true
 }
}


// print log info

log.info ""
log.info "Functional annotation pipeline"
log.info ""
log.info "General parameters"
log.info "------------------"

if ( params.proteinFile == null || params.proteinFile == "" ) {
  log.info "No protein sequence file specified!"
  exit 1
} else {
  if ( file( params.proteinFile ).exists() && file( params.proteinFile ).size() > 0 ) {
    log.info "Protein sequence file            : ${params.proteinFile}"
    protein = file(params.proteinFile)
  } else {
    log.info "Protein sequence file does not exist or it is empty!"
    exit 1
  }
}

if ( params.gffFile == null || params.gffFile == "" ) {
  log.info "No GFF Structural Annotation file specified!"
  log.info "We proceed anyway..."
} else {
  if ( file( params.gffFile ).exists() && file( params.gffFile ).size() > 0 ) {
    log.info "GFF Structural Annotation file              : ${params.gffFile}"
    gffavail = true
    annotation = file(params.gffFile)
  } else {
    log.info "GFF Structural Annotation file is missing or empty."
    log.info "We stop the pipeline so you can check it and define as \"\" otherwise if no GFF file is provided."
    exit 1
  }
}

if ( params.blastFile ) {
  log.info "BLAST results file           : ${params.blastFile}"
}

log.info "Species name                  : ${params.speciesName}"
log.info "KEGG species                  : ${params.kegg_species}"

if ( mysql ) {
  log.info "MySQL FA database 		       : ${params.dbname}"
} else {
  log.info "SQLite FA database 		       : $dbFileName"
}

if ( skip_cdSearch ) {
  log.info "CD-Search queries will be skipped."
}

if ( skip_sigtarp ) {
  log.info "SignalP and targetP queries will be skipped."
}

// split protein fasta file into chunks and then execute annotation for each chunk
// chanels for: interpro, blast, signalP, targetP, cdsearch_hit, cdsearch_features

chunkSize = params.chunkSize
chunkBlastSize = chunkSize
chunkIPSSize = chunkSize
chunkKoalaSize = chunkSize
chunkWebSize = chunkSize

if ( params.chunkBlastSize ) {
  chunkBlastSize = params.chunkBlastSize
}

if ( params.chunkIPSSize ) {
  chunkIPSSize = params.chunkIPSSize
}

if ( params.chunkKoalaSize ) {
  chunkKoalaSize = params.chunkKoalaSize
}

if ( params.chunkWebSize ) {
  chunkWebSize = params.chunkWebSize
}

workflow {
    FA-NF ()
}
