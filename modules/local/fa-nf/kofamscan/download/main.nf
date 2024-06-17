process KOFAMSCAN_DOWNLOAD {

    tag "$meta.id"
    label 'process_single'

    container 'docker.io/guigolab/fa-nf:0.4.0'

    input:
    tuple val(meta), val(ko_version)

    output:
    tuple val(meta), path("ko_list"),  emit: ko_list
    tuple val(meta), path("profiles"), emit: profiles
    tuple val(meta), path("ko_store"), emit: ko_store
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    ko_url = "https://www.genome.jp/ftp/db/kofam/archives"
    """
    curl -L --retry 3 -o ko_list.gz ${ko_url}/${ko_version}/ko_list.gz;
    gunzip ko_list.gz;
    curl -L --retry 3 -o profiles.tar.gz ${ko_url}/${ko_version}/profiles.tar.gz;
    tar zxf profiles.tar.gz; rm profiles.tar.gz;
    mkdir ko_store
    bulkDownloadKEGG.pl ko_list ko_store

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bulkDownloadKEGG: 0.0.0
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch ko_list
    mkdir profiles
    mkdir ko_store

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bulkDownloadKEGG: 0.0.0
    END_VERSIONS
    """

}
