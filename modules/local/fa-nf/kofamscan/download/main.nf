process KOFAMSCAN_DOWNLOAD {

    publishDir params.dbKOPath, mode: 'copy'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/interproscan:5.59_91.0--hec16e2b_1' :
        'biocontainers/interproscan:5.59_91.0--hec16e2b_1' }"

    label 'download'

    output:
    file "ko_list" into ko_list
    file "profiles" into ko_profiles
    file "ko_store" into ko_store

    """
    curl --retry 3 -o ko_list.gz ${params.koURLlist};
    gunzip ko_list.gz;
    curl --retry 3 -o profiles.tar.gz ${params.koURLprofiles};
    tar zxf profiles.tar.gz; rm profiles.tar.gz;
    mkdir ko_store
    bulkDownloadKEGG.pl ko_list ko_store
    """

}
