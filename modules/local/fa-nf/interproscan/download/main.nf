process INTERPROSCAN_DOWNLOAD {
    tag "$fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/interproscan:5.59_91.0--hec16e2b_1' :
        'biocontainers/interproscan:5.59_91.0--hec16e2b_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path ("data"), emit: data_interpro
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    curl --retry 3 -o iprscan.tar.gz ${params.iprscanURL};
    tar zxf iprscan.tar.gz
    rm iprscan.tar.gz
    cd interproscan-${params.iprscanVersion}
    python3 initial_setup.py
    cd ..
    mv interproscan-${params.iprscanVersion}/data .
    rm -rf interproscan-${params.iprscanVersion}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \$(echo \$(interproscan.sh --version 2>&1) | head -n1 | sed 's/^InterProScan version //;')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p data
    touch data/dummy
    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        interproscan: \$(echo \$(interproscan.sh --version 2>&1) | head -n1 | sed 's/^InterProScan version //;')
    END_VERSIONS
    """
}
