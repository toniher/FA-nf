process KOFAMSCAN_PARSE {
    tag "$fasta"
    label 'process_low'

    container ""

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("out_signalp_*"), emit: out_signalp
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    signalp -fasta ${fasta} \
    -stdout > out_signalp_${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        signalp: \$(echo \$(signalp -version --version 2>&1) | sed 's/^.*SignalP version //; s/Linux.*$//')
    END_VERSIONS
    """

    stub:
    """
    touch out_signalp_${fasta}
    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        signalp: \$(echo \$(signalp -version --version 2>&1) | sed 's/^.*SignalP version //; s/Linux.*$//')
    END_VERSIONS
    """
}
