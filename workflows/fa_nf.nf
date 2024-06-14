#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// TODO: Validate inputs here

//  Include modules below
include { AGAT_CONVERTSPGXF2GXF } from '../modules/nf-core/agat/convertspgxf2gxf/main'
include { DIAMOND_BLASTP } from '../modules/nf-core/diamond/blastp/main'
include { KOFAMSCAN } from '../modules/nf-core/kofamscan/main'
include { INTERPROSCAN } from '../modules/nf-core/interproscan/main'


mysql = false

workflow FA_NF {

}

// On finising
workflow.onComplete {

    println ( workflow.success ? "\nDone! Check results in --> $params.outdir\n" : "Oops .. something went wrong" )

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
        if (params.email) {
            log.info "Sending email to ${params.email}\n"
            sendMail(to: params.email, subject: "[FA-nf] Execution finished", body: msg)
        }
    }
}

workflow.onError {

    println( "Something went wrong" )

    if ( mysql ) {

    def procfile = new File( params.mysqllog+"/PROCESS" )
        procfile.delete()
    }

}
