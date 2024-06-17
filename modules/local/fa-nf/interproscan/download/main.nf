process INTERPROSCAN_DOWNLOAD {

    tag "$meta.id"
    label 'process_low'

    container 'docker.io/biocorecrg/interproscan:${interproscan_version}' // TODO: To fix

    input:
    tuple val(meta), val(interproscan_version)

    output:
    tuple val(meta), path("data"),     emit: interproscan_data
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    interproscan_url = "https://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/5/${params.interproscan_version}/alt/interproscan-data-${interproscan_version}.tar.gz"
    """
    curl -L --retry 3 -o iprscan.tar.gz ${interproscan_url};
    tar zxf iprscan.tar.gz
    rm iprscan.tar.gz
    cd interproscan-${params.interproscan_version}
    python3 initial_setup.py
    cd ..
    mv interproscan-${params.interproscan_version}/data .
    rm -rf interproscan-${params.interproscan_version}

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
