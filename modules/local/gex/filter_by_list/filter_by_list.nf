process FILTER_BY_LIST {
    label 'process_single'
    scratch true

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'biocontainers/pandas:1.5.2' }"

    input:
    val dataset_name
    path data_path
    path genes_path

    output:
    path "${dataset_name}.filtered_by_list.csv", emit: filtered_by_list_csv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def genes_file = genes_path ?: ''
    """
    #!/usr/bin/env python3

    import pandas as pd
    import sys

    def filter_by_list():
        data = pd.read_csv("${data_path}", index_col=0)
        genes_file = "${genes_file}"

        if genes_file:
            with open(genes_file) as handle:
                genes = [line.strip() for line in handle if line.strip()]

            data = data.loc[data.index.isin(genes)]

        data.to_csv("${dataset_name}.filtered_by_list.csv")

    filter_by_list()

    # Create versions file
    import pandas
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f'    python: "{sys.version.split()[0]}"\\n')
        f.write(f'    pandas: "{pandas.__version__}"\\n')
    """
}
