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

include { FA_NF } from './workflows/fa_nf'
include { DOWNLOAD } from './workflows/download'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_fanf_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_fanf_pipeline'

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

// print log info

log.info ""
log.info "Functional annotation pipeline"
log.info ""
log.info "General parameters"
log.info "------------------"

workflow RUN_DOWNLOAD {

    take:
    dbnames // channel: dbnames from --input

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    DOWNLOAD(dbnames)
}

workflow RUN_FA_NF {

    main:
    FA_NF (params.fasta)
}

workflow {

    main:
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.dblist
    )

    RUN_DOWNLOAD (
        PIPELINE_INITIALISATION.out.dbnames
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    // PIPELINE_COMPLETION (
    //     params.email,
    //     params.email_on_fail,
    //     params.plaintext_email,
    //     params.outdir,
    //     params.monochrome_logs,
    //     params.hook_url
    // )

}
