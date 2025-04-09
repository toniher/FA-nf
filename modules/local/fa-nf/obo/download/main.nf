process OBO_DOWNLOAD {

    tag "$meta.id"
    label 'process_single'

    container 'debian-perlbrew:latest'

    input:
    tuple val(meta), val(url)

    output:
    tuple val(meta), path("*.obo"),  emit: obo_file

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: meta.id
    """
    curl -L --retry 3 -o gene_ontology.obo ${url};

    """

    stub:
    prefix = task.ext.prefix ?: meta.id
    """
    touch gene_ontology.obo

    """

}
