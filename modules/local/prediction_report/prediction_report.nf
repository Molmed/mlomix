process PREDICTION_REPORT {
    label 'process_single'
    publishDir "${params.outdir}/finalize", mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    conda "${moduleDir}/environment.yml"

    input:
    path predictions
    path labels
    path colors
    val modality

    output:
    path "prediction_sets.${modality}.upset.png", emit: upset_plot
    path "prediction_sets.${modality}.sankey.html", emit: sankey_plot
    path "versions.yml", emit: versions

    script:
    """
    python "${projectDir}/bin/prediction_plots.py" \
        --predictions "${predictions}" \
        --labels "${labels}" \
        --colors "${colors}" \
        --modality "${modality}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.11"
        pandas: "2.2.3"
        matplotlib: "3.9.2"
        upsetplot: "0.9.0"
        plotly: "6.8.0"
    END_VERSIONS
    """
}
