process INFERENCE {
    label 'process_single'
    publishDir "${params.outdir}/finalize", mode: params.publish_dir_mode,
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }

    conda "${moduleDir}/environment.yml"

    input:
    path transposed
    path model
    path calibration
    val modality
    val alpha

    output:
    path "predictions.${modality}.csv", emit: predictions
    path "versions.yml", emit: versions

    script:
    """
    python "${projectDir}/bin/inference.py" \
        --input "${transposed}" \
        --model "${model}" \
        --calibration "${calibration}" \
        --output "predictions.${modality}.csv" \
        --alpha "${alpha}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "3.11"
        joblib: "1.4.2"
        numpy: "1.26.4"
        pandas: "2.2.3"
        scikit-learn: "1.5.2"
    END_VERSIONS
    """
}
