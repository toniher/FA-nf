#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// TODO: Validate inputs here

//  Include modules below
include { AGAT_CONVERTSPGXF2GXF } from './modules/nf-core/agat/convertspgxf2gxf/main'
include { DIAMOND_BLASTP } from '../modules/nf-core/diamond/blastp/main' 
include { KOFAMSCAN } from './modules/nf-core/kofamscan/main'
include { INTERPROSCAN } from '../modules/nf-core/interproscan/main'  

worflow FA-NF {

    seqData = Channel
     .from(protein)
     .splitFasta( by: chunkSize )

     seqBlastData = Channel
      .from(protein)
      .splitFasta( by: chunkBlastSize )

    seqKoalaData = Channel
     .from(protein)
     .splitFasta( by: chunkKoalaSize )

     seqIPSData = Channel
      .from(protein)
      .splitFasta( by: chunkIPSSize )

    seqWebData = Channel
     .from(protein)
     .splitFasta( by: chunkWebSize )

    ipscan_properties = file(params.ipscanproperties)

    if ( params.debug == "true" || params.debug == true ) {
     println("Debugging... only the first $params.debugSize chunks will be processed")
     // Diferent parts for different processes.
     // TODO: With DSL2 this is far simpler
     (seq_file1, seq_file2) = seqData.take(params.debugSize).into(2)
     (seq_file_blast) = seqBlastData.take(params.debugSize).into(1)
     (seq_file_koala) = seqKoalaData.take(params.debugSize).into(1)
     (seq_file_ipscan) = seqIPSData.take(params.debugSize).into(1)
     (web_seq_file1, web_seq_file2) = seqWebData.take(params.debugSize).into(2)

     testNum = ( params.chunkSize.toInteger() * params.debugSize )
     seqTestData = Channel
      .from(protein)
      .splitFasta(by: testNum)

      (seq_test) = seqTestData.take(1).into(1)

    } else {
     println("Process entire dataset")
     (seq_file1, seq_file2) = seqData.into(2)
     (seq_file_blast) = seqBlastData.into(1)
     (seq_file_koala) = seqKoalaData.into(1)
     (seq_file_ipscan) = seqIPSData.into(1)
     (web_seq_file1, web_seq_file2) = seqWebData.into(2)

     seqTestData = Channel
      .from(protein)

     // Anything for keeping. This is only kept for coherence
     (seq_test) = seqTestData.into(1)

    }

    if ( params.oboFile == "" || params.oboFile == null ) {
      oboFile = downloadURL( "http://www.geneontology.org/ontology/gene_ontology.obo", "gene_ontology.obo" )
    } else {
      oboFile = params.oboFile
    }


    // Preprocessing GFF File
    if ( gffavail ) {

      if ( gffclean ) {

       process cleanGFF {

        publishDir params.resultPath, mode: 'copy'

        label 'gffcheck'

        input:
         file config_file

        output:
         file "annot.gff" into gff_file
         file "annot.gff.clean.txt" into gff_file_log

         """
          # get annot file
          export escaped=\$(echo '$baseDir')
          export basedirvar=\$(echo '\\\$\\{baseDir\\}')
          export basedirvar=\$(echo '\\\$\\{baseDir\\}')
          input_file=`perl -lae 'if (\$_=~/gffFile\\s*\\=\\s*[\\x27|\\"](\\S+)[\\x27|\\"]/) { \$base = \$1; \$base=~s/\$ENV{'basedirvar'}/\$ENV{'escaped'}/g; print \$base }' $config_file`

          if [ "\${input_file: -3}" == ".gz" ]; then
            gunzip -c \$input_file > input_gff
            input_file=input_gff
          fi
          agat_convert_sp_gxf2gxf.pl --gff \$input_file -o annot.gff > annot.gff.clean.txt
         """

       }


      } else {

       process copyGFF {

        label 'gffcheck'

        input:
         file config_file

        output:
         file "annot.gff" into gff_file

         """
          # get annot file
          export escaped=\$(echo '$baseDir')
          export basedirvar=\$(echo '\\\$\\{baseDir\\}')
          input_file=`perl -lae 'if (\$_=~/gffFile\\s*\\=\\s*[\\x27|\\"](\\S+)[\\x27|\\"]/) { \$base = \$1; \$base=~s/\$ENV{'basedirvar'}/\$ENV{'escaped'}/g; print \$base }' $config_file`

          if [ "\${input_file: -3}" == ".gz" ]; then
            gunzip -c \$input_file > input_gff
            input_file=input_gff
          fi

          cp \$input_file annot.gff
         """

       }
      }

      if ( gffstats ) {

       process statsGFF {

        publishDir params.resultPath, mode: 'copy'

        label 'gffcheck'

        input:
         file gff_file

        output:
         file "*.txt" into gff_stats

         """
          # Generate Stats
          agat_sp_statistics.pl --gff $gff_file > ${gff_file}.stats.txt
         """

       }


      }

    } else {

      // Dummy empty GFF
      process dummyGFF {

       label 'gffcheck'

       input:
        file config_file

       output:
        file "annot.gff" into gff_file

        """
         # empty annot file
         touch annot.gff
        """

      }
    }

    // Blast like processes
    // TODO: To change for different aligners
    diamond = false

    if( params.diamond == "true" || params.diamond == true ) {
     diamond = true
    }

    // BlastAnnotMode
    blastAnnotMode = "common"
    if( params.blastAnnotMode != "" && params.blastAnnotMode != null ) {
      blastAnnotMode = params.blastAnnotMode
    }

    if ( params.blastFile == "" ||  params.blastFile == null ){

     // program-specific parameters
     db_name = file(blastDbPath).name
     db_path = file(blastDbPath).parent

     // Handling Database formatting
     formatdbDetect = "false"

     if ( diamond ) {

      formatDbFileName = blastDbPath + ".dmnd"
      formatDbFile = file(formatDbFileName)
      if ( formatDbFile.exists() && formatDbFile.size() > 0 ) {
       formatdbDetect = "true"
      }

      if ( formatdbDetect == "false" ) {

       process diamondFormat{

        label 'diamond'

        output:
        file "${db_name}_formatdb.dmnd" into formatdb

        """
         diamond makedb --in ${db_path}/${db_name} --db "${db_name}_formatdb"
        """
       }

      } else {
       formatdb = blastDbPath
      }

     } else {

      formatDbDir = file( db_path )
      filter =  ~/${db_name}.*.phr/
      def fcount = 0
      formatDbDir.eachFileMatch( filter ) { it ->
       fcount = fcount + 1
      }
      if ( fcount > 0 ) {
        formatdbDetect = "true"
      }

      println( formatdbDetect )
      if ( formatdbDetect == "false" ) {

       // println( "TUR" )

       process blastFormat{

        label 'blast'

        output:
        file "${db_name}.p*" into formatdb

        """
         makeblastdb -dbtype prot -in ${db_path}/${db_name} -parse_seqids -out ${db_name}
        """
       }

      } else {
       formatdb = blastDbPath
      }
     }

     if ( diamond == true ) {

      process diamond{

       label 'diamond'

       input:
       file seq from seq_file_blast
       file formatdb_file from formatdb

       output:
       file "blastXml${seq}" into (blastXmlResults1, blastXmlResults2, blastXmlResults3)

       script:
       if ( formatdbDetect == "false" ) {
        command = "diamond blastp --db ${formatdb_file} --query $seq --outfmt 5 --threads ${task.cpus} --evalue ${params.evalue} --out blastXml${seq}"
       } else {
        command = "diamond blastp --db ${db_path}/${db_name} --query $seq --outfmt 5 --threads ${task.cpus} --evalue ${params.evalue} --out blastXml${seq}"
       }

       command

      }

     } else {

      process blast{

       label 'blast'

       // publishDir "results", mode: 'copy'

       input:
       file seq from seq_file_blast
       file formatdb_file from formatdb

       output:
       file "blastXml${seq}" into (blastXmlResults1, blastXmlResults2, blastXmlResults3)

       script:
       if ( formatdbDetect == "false" ) {
        command = "blastp -db ${formatdb_file} -query $seq -num_threads ${task.cpus} -evalue ${params.evalue} -out blastXml${seq} -outfmt 5"
       } else {
        command = "blastp -db ${db_path}/${db_name} -query $seq -num_threads ${task.cpus} -evalue ${params.evalue} -out blastXml${seq} -outfmt 5"
       }

       command
      }

     }

    } else {

     blastInput=file(params.blastFile)

     process convertBlast {

      // publishDir "results", mode: 'copy'

      input:
      file blastFile from blastInput

      output:
      file("*.xml") into (blastXmlResults1, blastXmlResults2, blastXmlResults3)

      """
       hugeBlast2XML.pl -blast $blastFile -n 1000 -out blast.res
      """

     }
    }

    if ( kolist != "" ||  kolist != null ){

      process kofamscan{

       label 'kofamscan'

       input:
       file seq from seq_file_koala

       output:
       file "koala_${seq}" into koalaResults

       """
        exec_annotation --cpu ${task.cpus} -p ${koprofiles} -k ${kolist} -o koala_${seq} $seq
       """

      }

      process kofam_parse {

       input:
       file "koala_*" from koalaResults.collect()

       output:
       file allKoala into ( koala_parsed, koala_parsed2 )

      """

      mkdir -p output
      processHmmscan2TSV.pl "koala_*" output
      cat output/koala_* > allKoala
      """

      }

      // Replacing keggfile
      keggfile = koala_parsed

    } else {


     if (params.keggFile == "" ||  params.keggFile == null ) {

      println "Please run KEGG KO group annotation on the web server http://www.genome.jp/tools/kaas/"
      exit 1
     } else {

       keggfile = file(params.keggFile)
     }

    }

    // GO retrieval from BLAST results
    if (params.gogourl != "") {

      process blast_annotator {

       label 'blastannotator'

       maxForks 3

       input:
       file blastXml from blastXmlResults2.flatMap()

       output:
       file "blastAnnot" into ( blast_annotator_results1, blast_annotator_results2 )

      """
       blast-annotator.pl -in $blastXml -out blastAnnot --hits $params.gogohits --url $params.gogourl -t $blastAnnotMode -q --format blastxml
      """
      }

    }

    process blastDef {

     // publishDir "results", mode: 'copy'
     tag "${blastXml}"

     input:
     file blastXml from blastXmlResults3.flatMap()

     output:
     file "blastDef_${blastXml}.txt" into blastDef_results

     """
      definitionFromBlast.pl  -in $blastXml -out blastDef_${blastXml}.txt -format xml -q
     """
    }

    process ipscn {

        label 'ipscan'

        if ( workflow.containerEngine == "singularity" ) {
          containerOptions "--bind ${ipscandata}:/usr/local/interproscan/data"
        } else {
          containerOptions "--volume ${ipscandata}:/usr/local/interproscan/data"
        }

        input:
        file seq from seq_file_ipscan
        file ("interproscan.properties") from file( ipscan_properties )

        output:
        file("out_interpro_${seq}") into (ipscn_result1, ipscn_result2)

        """
        sed 's/*//g' $seq > tmp4ipscn
        interproscan.sh -i tmp4ipscn --goterms --iprlookup --pathways -o out_interpro_${seq} -f TSV -T ${params.ipscantmp}
        """
    }

    process 'cdSearchHit' {

        label 'cdSearch'

        if ( ! skip_cdSearch ) {
          maxForks 1
        }

        input:
        file seq from web_seq_file1

        output:
        file("out_hit_${seq}") into ( cdSearch_hit_result1, cdSearch_hit_result2 )

        script:
        if ( skip_cdSearch ) {
          // Dummy content
          command = "touch out_hit_${seq}"
        } else {
          command = "submitCDsearch.pl -o out_hit_${seq} -in $seq"
        }

        command
    }

    process 'cdSearchFeat' {

        label 'cdSearch'

        if ( ! skip_cdSearch ) {
          maxForks 1
        }

        input:
        file seq from web_seq_file2

        output:
        file("out_feat_${seq}") into ( cdSearch_feat_result1, cdSearch_feat_result2 )

        script:
        if ( skip_cdSearch ) {
          // Dummy content
          command = "touch out_feat_${seq}"
        } else {
          command = "submitCDsearch.pl -t feats -o out_feat_${seq} -in $seq"
        }

        command
    }

    if ( skip_sigtarp ) {

      process 'signalP_dummy' {

          input:
          file seq from seq_file1

          output:
          file("out_signalp_${seq}") into (signalP_result1, signalP_result2)

          """
          touch out_signalp_${seq}
          """
      }

      process 'targetP_dummy' {

          input:
          file seq from seq_file2

          output:
          file("out_targetp_${seq}") into (targetP_result1, targetP_result2)

          """
          touch out_targetp_${seq}
          """
      }

    } else {


      process 'signalP' {

          label 'sigtarp'

          input:
          file seq from seq_file1

          output:
          file("out_signalp_${seq}") into (signalP_result1, signalP_result2)

          """
          signalp -fasta $seq -stdout > out_signalp_${seq}
          """
      }

      process 'targetP' {

          label 'sigtarp'

          input:
          file seq from seq_file2

          output:
          file("out_targetp_${seq}") into (targetP_result1, targetP_result2)

          """
          targetp -fasta $seq -stdout > out_targetp_${seq}
          """
      }

    }

    // Database setup below
    process initDB {

     input:
      file config_file
      file gff_file
      file seq from seq_test

     output:
      file 'config' into (config4perl7, config4perl8, config4perl10)

     script:
     command = "mkdir -p $params.resultPath\n"
     command += "sed 's/^\\s*params\\s*{\\s*\$//gi' $config_file | sed 's/^\\s*}\\s*\$//gi' | sed '/^\\s*\$/d' | sed 's/\\s\\=\\s/:/gi' | sed '/^\\s*\\/\\//d' > configt\n"
     command += "export escaped=\$(echo '$baseDir')\n"
     command += "export basedirvar=\$(echo '\\\$\\{baseDir\\}')\n"
     command += "perl -lae '\$_=~s/\$ENV{'basedirvar'}/\$ENV{'escaped'}/g; print;' configt > config\n"


     if ( mysql ) {
      // Add dbhost to config
      command += "echo \"\$(cat config)\n dbhost:${dbhost}\" > configIn ;\n"
      command += "fa_main.v1.pl init -conf configIn"

       if ( gffavail && gffclean ) {
        command += " -gff ${gff_file}"
       }
     } else {

        if (exists) {
         log.info "SQLite database ${dbFileName} exists. We proceed anyway..."
        }

        command += "fa_main.v1.pl init -conf config"

        if ( gffavail && gffclean ) {
         command += " -gff ${gff_file}"
        }
     }

     if ( params.debug=="TRUE"||params.debug=="true" ) {
       // If in debug mode, we restrict de seq entries we process
       command += " -fasta ${seq}"
     }

     if ( params.rmversion=="TRUE"||params.rmversion=="true" ) {
       // If remove versioning in protein sequences (for cases like ENSEMBL)
       command += " -rmversion"
     }

     command
    }


    if ( ! koentries ) {

      process 'kegg_download'{

       maxForks 1

       input:
       file keggfile from keggfile
       file config from config4perl8

       output:
       file("down_kegg") into (down_kegg)


       script:

        command = "download_kegg_KAAS.pl -input $keggfile -conf $config > done 2>err"

        command
      }

    } else {

      process 'kegg_download_dummy' {

       input:
       file keggfile from keggfile
       file config from config4perl8

       output:
       file("down_kegg") into (down_kegg)


       script:

        command = "touch down_kegg"

        command

      }
    }

    // Data upload process
    process 'data_upload' {

     maxForks 1

     label 'upload'

     input:

      file "def*" from blastDef_results.collect()

      file "out_signalp*" from signalP_result1.collect()
      file "out_targetp*" from targetP_result1.collect()

      file "out_interpro*" from ipscn_result1.collect()

      file "out_hit*" from cdSearch_hit_result1.collect()
      file "out_feat*" from cdSearch_feat_result1.collect()

      file keggfile from keggfile

      file("down_kegg") from down_kegg

      file "blastAnnot*" from blast_annotator_results1.collect()

      file config from config4perl7

      output:
      file('done') into (last_step)

     script:

      command = checkMySQL( mysql, params.mysqllog )

      command += " \
       cat def* > allDef; \
       upload_go_definitions.pl -i allDef -conf \$config -mode def -param 'blast_def' > def_done ; \
      "

      command += " \
       cat out_signalp* > allSignal ; \
       load_sigtarp.pl -i allSignal -conf \$config -type s > upload_signalp ; \
      "

      command += " \
       cat out_targetp* > allTarget ; \
       load_sigtarp.pl -i allTarget -conf \$config -type t > upload_targetp ; \
      "

      command += " \
       cat out_interpro* > allInterpro ; \
       run_interpro.pl -mode upload -i allInterpro -conf \$config > upload_interpro ; \
      "

      command += " \
       cat out_hit* > allCDsearchHit ; \
       upload_CDsearch.pl -i allCDsearchHit -type h -conf \$config > upload_hit ; \
      "

      command += " \
       cat out_feat* > allCDsearchFeat ; \
       upload_CDsearch.pl -i allCDsearchFeat -type f -conf \$config > upload_feat ; \
      "

      command += " \
       cat blastAnnot* > allBlast ; \
       awk '\$2!=\"#\"{print \$1\"\t\"\$2}' allBlast > two_column_file_blast ; \
       upload_go_definitions.pl -i two_column_file_blast -conf \$config -mode go -param 'blast_annotator' > done ; \
      "

      if ( ! koentries ) {
        command += " \
         load_kegg_KAAS.pl -input $keggfile -dir down_kegg -rel $params.kegg_release -conf \$config > upload_kegg 2>err; \
        "
      } else {
        command += " \
         load_kegg_KAAS.pl -input $keggfile -entries $koentries -rel $params.kegg_release -conf \$config > upload_kegg 2>err; \
        "
      }

      command
    }

    process 'generateResultFiles'{
     input:
      file config from config4perl10
      file all_done from last_step

     script:

      command = checkMySQL( mysql, params.mysqllog )

      if ( annotation != null && annotation != "" ){
        command += " \
         get_gff3.pl -conf \$config ; \
        "
      }

      command += " \
       get_results.pl -conf \$config -obo ${oboFile} ; \
      "

      command
    }

    // Check MySQL IP
    def checkMySQL( mysql, mysqllog )  {

     command = ""

     if ( mysql ) {
       // Add dbhost to config
       command += "DBHOST=\"dbhost:'`cat ${mysqllog}/DBHOST`'\"; echo \"\$(cat config)\n \$DBHOST\" > configIn ;\n"
       command += "config=configIn ;"
     } else {
       command += "config=config ;"
     }

     return command

    }

    def downloadURL( address, filename ) {
      downFile = new File( filename ) << new URL (address).getText()
      return downFile.absolutePath
    }


    if ( ! skip_sigtarp ) {
      signalP_result2
       .collectFile(name: file(params.resultPath + "signalP.res.tsv"))
        .println { "Result saved to file: $it" }

      targetP_result2
       .collectFile(name: file(params.resultPath + "targetP.res.tsv"))
        .println { "Result saved to file: $it" }
    }

    if ( ! skip_cdSearch ) {

      cdSearch_hit_result2
       .collectFile(name: file(params.resultPath + "cdSearch_hit.res.tsv"))
        .println { "Result saved to file: $it" }

      cdSearch_feat_result2
       .collectFile(name: file(params.resultPath + "cdSearch_feat.res.tsv"))
        .println { "Result saved to file: $it" }

    }

    ipscn_result2
      .collectFile(name: file(params.resultPath + "interProScan.res.tsv"))
      .println { "Result saved to file: $it" }

    if ( kolist != "" ||  kolist != null ){

      koala_parsed2
        .collectFile(name: file(params.resultPath + "koala.res.tsv"))
        .println { "Result saved to file: $it" }
    }

    if (params.gogourl != "") {
      blast_annotator_results2
        .collectFile(name: file(params.resultPath + "blastAnnotator.res.tsv"))
        .println { "Result saved to file: $it" }
    }

}
// On finising
workflow.onComplete {

    println ( workflow.success ? "\nDone! Check results in --> $params.resultPath\n" : "Oops .. something went wrong" )

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

    if ( mysql ) {

    def procfile = new File( params.mysqllog+"/PROCESS" )
        procfile.delete()
    }

}
