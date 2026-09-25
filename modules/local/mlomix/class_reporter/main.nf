process CLASS_FILTER_AND_REPORT {
    cache 'lenient'
    tag "class_filter_and_report"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "biocontainers/pandas:1.5.3_cv1"

    input:
    path samplesheet
    val gex_classes_file
    val dnam_classes_file

    output:
    path "class_distribution.gex.all.svg"     , emit: class_plot_gex_all
    path "class_distribution.gex.filtered.svg", emit: class_plot_gex_filtered
    path "class_distribution.dnam.all.svg"    , emit: class_plot_dnam_all
    path "class_distribution.dnam.filtered.svg", emit: class_plot_dnam_filtered
    path "classes.gex.all.tsv"                , emit: classes_gex_all
    path "classes.gex.filtered.tsv"           , emit: classes_gex_filtered
    path "classes.dnam.all.tsv"               , emit: classes_dnam_all
    path "classes.dnam.filtered.tsv"          , emit: classes_dnam_filtered
    path "samplesheet.filtered.csv"           , emit: samplesheet_filtered_csv
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gex_class_file_arg = gex_classes_file ? "--gex-classes-file ${gex_classes_file}" : ""
    def dnam_classes_file_arg = dnam_classes_file ? "--dnam-classes-file ${dnam_classes_file}" : ""
    def gex_min_arg = params.gex_min_samples_per_class != null ? "--gex-min ${params.gex_min_samples_per_class}" : ""
    def dnam_min_arg = params.dnam_min_samples_per_class != null ? "--dnam-min ${params.dnam_min_samples_per_class}" : ""
    """
    python3 ${projectDir}/bin/class_filter_and_report.py ${samplesheet} ${gex_min_arg} ${dnam_min_arg} ${gex_class_file_arg} ${dnam_classes_file_arg}
    """

    stub:
    """
    touch class_distribution.gex.all.svg
    touch class_distribution.gex.filtered.svg
    touch class_distribution.dnam.all.svg
    touch class_distribution.dnam.filtered.svg
    touch classes.gex.all.tsv
    touch classes.gex.filtered.tsv
    touch classes.dnam.all.tsv
    touch classes.dnam.filtered.tsv
    touch samplesheet.filtered.csv
    echo "versions:" > versions.yml
    """
}
