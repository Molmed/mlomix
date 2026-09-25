/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GEX_REF_PREPROCESSOR } from './gex_ref_preprocessor'
include { GEX } from './gex'
include { DNAM } from './dnam'
include { FINALIZE } from '../modules/local/finalize/finalize'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MLOMIX {

    take:
    ch_gex_samplesheet
    ch_datasets
    ch_batches
    ch_classes
    ch_class_plot_gex_filtered
    ch_class_plot_dnam_filtered
    ch_classes_gex_filtered
    ch_classes_dnam_filtered
    ch_samplesheet_filtered_csv
    ch_dnam_samplesheet
    ch_dnam_beta_matrix
    ch_dnam_pvals
    classifier_config
    genome
    annotation_version
    random_seed
    run_gex
    run_dnam
    use_precomputed_dnam

    main:

    ch_versions = channel.empty()
    ch_gex_transposed = channel.empty()
    ch_gex_norm_factors = channel.empty()
    ch_dnam_transposed = channel.empty()
    ch_classes_tsv = channel.empty()
    ch_visuals = channel.empty()

    ch_genome_for_gex = genome
        .combine(run_gex)
        .filter { item -> item[1] }
        .map { item -> item[0] }

    ch_annotation_for_gex = annotation_version
        .combine(run_gex)
        .filter { item -> item[1] }
        .map { item -> item[0] }

    ch_gex_samplesheet_gated = ch_gex_samplesheet
        .combine(run_gex)
        .filter { item -> item[-1] }
        .map { item -> [item[0], item[1]] }

    ch_datasets_gated = ch_datasets
        .combine(run_gex)
        .filter { item -> item[-1] }
        .map { item -> [item[0], item[1], item[2]] }

    ch_batches_gated = ch_batches
        .combine(run_gex)
        .filter { item -> item[-1] }
        .map { item -> item[0] }

    ch_classes_gated = ch_classes
        .combine(run_gex)
        .filter { item -> item[-1] }
        .map { item -> item[0] }

    ch_random_seed_gated = run_gex
        .filter { _flag -> _flag }
        .map { _flag -> random_seed }

    GEX_REF_PREPROCESSOR (
        ch_genome_for_gex,
        ch_annotation_for_gex
    )

    GEX (
        ch_gex_samplesheet_gated,
        ch_datasets_gated,
        ch_batches_gated,
        ch_classes_gated,
        GEX_REF_PREPROCESSOR.out.filtered_annotations,
        classifier_config.gex_genes_file,
        classifier_config.gex_norm_factors_file ?: '',
        classifier_config.class_colors_file ?: '',
        ch_random_seed_gated
    )

    ch_versions = ch_versions.mix(GEX_REF_PREPROCESSOR.out.versions)
    ch_versions = ch_versions.mix(GEX.out.versions)
    ch_gex_transposed = GEX.out.transposed_csv
    ch_gex_norm_factors = GEX.out.norm_factors_rds
    ch_classes_tsv = ch_classes
    ch_visuals = ch_visuals.mix(GEX.out.visuals)

    // Pass ch_dnam_samplesheet directly without restructuring to preserve list integrity
    ch_dnam_samplesheet_gated = ch_dnam_samplesheet

    ch_use_precomputed_dnam_gated = use_precomputed_dnam
        .combine(run_dnam)
        .filter { item -> item[-1] }
        .map { item -> item[0] }

    DNAM (
        ch_dnam_samplesheet_gated,
        ch_classes,
        random_seed,
        ch_dnam_beta_matrix,
        ch_dnam_pvals,
        ch_use_precomputed_dnam_gated,
        classifier_config.dnam_probes_file,
        classifier_config.class_colors_file ?: ''
    )
    ch_versions = ch_versions.mix(DNAM.out.versions)
    ch_dnam_transposed = DNAM.out.transposed_csv.ifEmpty(file("${params.outdir}/.dnam_placeholder"))
    ch_gex_transposed = ch_gex_transposed.ifEmpty(file("${params.outdir}/.gex_placeholder"))
    ch_gex_norm_factors = ch_gex_norm_factors.ifEmpty(file("${params.outdir}/.gex_norm_factors_placeholder"))
    ch_visuals = ch_visuals.mix(DNAM.out.visuals)

    FINALIZE (
        ch_classes_tsv,
        ch_dnam_transposed,
        ch_gex_transposed,
        ch_gex_norm_factors,
        ch_visuals.toList(),
        ch_class_plot_gex_filtered,
        ch_class_plot_dnam_filtered,
        ch_classes_gex_filtered,
        ch_classes_dnam_filtered,
        ch_samplesheet_filtered_csv
    )
    ch_versions = ch_versions.mix(FINALIZE.out.versions)

    emit:
    versions = ch_versions
    labels_csv = FINALIZE.out.labels_csv
    dnam_features_csv = FINALIZE.out.dnam_features_csv
    gex_features_csv = FINALIZE.out.gex_features_csv
    gex_norm_factors_rds = FINALIZE.out.gex_norm_factors_rds
    classes_gex_filtered = FINALIZE.out.classes_gex_filtered
    classes_dnam_filtered = FINALIZE.out.classes_dnam_filtered
    samplesheet_filtered_csv = FINALIZE.out.samplesheet_filtered_csv
    class_distribution_gex_filtered_svg = FINALIZE.out.class_distribution_gex_filtered_svg
    class_distribution_dnam_filtered_svg = FINALIZE.out.class_distribution_dnam_filtered_svg
    visuals_dir = FINALIZE.out.visuals_dir

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
