#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/mlomix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/mlomix
    Website: https://nf-co.re/mlomix
    Slack  : https://nfcore.slack.com/channels/mlomix
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MLOMIX  } from './workflows/mlomix'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_mlomix_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_mlomix_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_mlomix_pipeline'
include { resolveClassifierConfig } from './subworkflows/local/utils_nfcore_mlomix_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Fetch the FASTA path from the selected iGenomes configuration.
params.fasta = getGenomeAttribute('fasta')
params.tsne = params.tsne ?: true
params.gex_norm_factors_file = params.gex_norm_factors_file ?: null

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_MLOMIX {

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

    //
    // WORKFLOW: Run pipeline
    //
    MLOMIX (
        ch_gex_samplesheet,
        ch_datasets,
        ch_batches,
        ch_classes,
        ch_class_plot_gex_filtered,
        ch_class_plot_dnam_filtered,
        ch_classes_gex_filtered,
        ch_classes_dnam_filtered,
        ch_samplesheet_filtered_csv,
        ch_dnam_samplesheet,
        ch_dnam_beta_matrix,
        ch_dnam_pvals,
        classifier_config,
        genome,
        annotation_version,
        random_seed,
        run_gex,
        run_dnam,
        use_precomputed_dnam
    )
    emit:
    versions = MLOMIX.out.versions
    classes_gex_filtered = MLOMIX.out.classes_gex_filtered
    classes_dnam_filtered = MLOMIX.out.classes_dnam_filtered
    samplesheet_filtered_csv = MLOMIX.out.samplesheet_filtered_csv
    class_distribution_gex_filtered_svg = MLOMIX.out.class_distribution_gex_filtered_svg
    class_distribution_dnam_filtered_svg = MLOMIX.out.class_distribution_dnam_filtered_svg
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    classifier_config = resolveClassifierConfig()
    classifier_config.each { key, value ->
        if (value != null) {
            params.put(key, value)
        }
    }

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
        classifier_config
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_MLOMIX (
        PIPELINE_INITIALISATION.out.gex_samplesheet,
        PIPELINE_INITIALISATION.out.datasets,
        PIPELINE_INITIALISATION.out.batches,
        PIPELINE_INITIALISATION.out.classes,
        PIPELINE_INITIALISATION.out.class_plot_gex_filtered,
        PIPELINE_INITIALISATION.out.class_plot_dnam_filtered,
        PIPELINE_INITIALISATION.out.classes_gex_filtered,
        PIPELINE_INITIALISATION.out.classes_dnam_filtered,
        PIPELINE_INITIALISATION.out.samplesheet_filtered_csv,
        PIPELINE_INITIALISATION.out.dnam_samplesheet,
        PIPELINE_INITIALISATION.out.dnam_beta_matrix,
        PIPELINE_INITIALISATION.out.dnam_pvals,
        classifier_config,
        PIPELINE_INITIALISATION.out.genome,
        PIPELINE_INITIALISATION.out.annotation_version,
        params.random_seed,
        PIPELINE_INITIALISATION.out.run_gex,
        PIPELINE_INITIALISATION.out.run_dnam,
        PIPELINE_INITIALISATION.out.use_precomputed_dnam
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        channel.empty()
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
