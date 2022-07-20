process TARGETP {
    tag "$fasta"
    label 'process_low'

    container ""

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("out_targetp_*"), emit: out_targetp
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    targetp -fasta ${fasta} \
    -stdout > out_targetp_${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        targetp: \$(echo \$(targetp -version --version 2>&1) | sed 's/^.*TargetP version //; s/Linux.*$//')
    END_VERSIONS
    """

    stub:
    """
    touch out_targetp_${fasta}
    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        targetp: \$(echo \$(targetp -version --version 2>&1) | sed 's/^.*TargetP version //; s/Linux.*$//')
    END_VERSIONS
    """
}
