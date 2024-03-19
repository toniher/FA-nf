process INTERPROSCAN_RUN {
    tag "$fasta"
    label 'process_low'

    container ""

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("out_interpro_*"), emit: out_interpro
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    sed 's/*//g' ${fasta} > tmp4ipscn

    interproscan.sh -i tmp4ipscn --goterms --iprlookup --pathways -o out_interpro_${fasta} -f TSV -T ${params.ipscantmp}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \$(echo \$(interproscan.sh --version 2>&1) | head -n1 | sed 's/^InterProScan version //;')
    END_VERSIONS
    """

    stub:
    """
    touch out_interpro_${fasta}
    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        interproscan: \$(echo \$(interproscan.sh --version 2>&1) | head -n1 | sed 's/^InterProScan version //;')
    END_VERSIONS
    """
}
