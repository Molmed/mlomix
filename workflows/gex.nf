/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { CONCATENATE_GEX        } from '../modules/local/gex/concatenate_gex'
include { UMAP as UMAP_GEX_RAW_BY_BATCH  } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_BC_BY_BATCH   } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_NORM_BY_BATCH } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_RAW_BY_CLASS  } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_BC_BY_CLASS   } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_NORM_BY_CLASS } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_PROCESSED_BY_BATCH } from '../modules/local/umap/umap'
include { UMAP as UMAP_GEX_PROCESSED_BY_CLASS } from '../modules/local/umap/umap'
include { TSNE as TSNE_GEX_RAW_BY_BATCH  } from '../modules/local/tsne/tsne'
include { TSNE as TSNE_GEX_BC_BY_BATCH   } from '../modules/local/tsne/tsne'
include { TSNE as TSNE_GEX_NORM_BY_BATCH } from '../modules/local/tsne/tsne'
include { TSNE as TSNE_GEX_RAW_BY_CLASS  } from '../modules/local/tsne/tsne'
include { TSNE as TSNE_GEX_BC_BY_CLASS   } from '../modules/local/tsne/tsne'
include { TSNE as TSNE_GEX_NORM_BY_CLASS } from '../modules/local/tsne/tsne'
include { TSNE as TSNE_GEX_PROCESSED_BY_BATCH } from '../modules/local/tsne/tsne'
include { TSNE as TSNE_GEX_PROCESSED_BY_CLASS } from '../modules/local/tsne/tsne'
include { BATCH_CORRECT          } from '../modules/local/gex/batch_correct'
include { FILTER_BY_ANNOTATION    } from '../modules/local/gex/filter_by_annotation/main'
include { FILTER_BY_LIST          } from '../modules/local/gex/filter_by_list/filter_by_list'
include { MERGE_DATASETS         } from '../modules/local/gex/merge_datasets'
include { NORMALIZE              } from '../modules/local/gex/normalize'
include { TRANSPOSE              } from '../modules/local/transpose'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow GEX {
    take:
    _ch_samplesheet       // channel: samplesheet from the samplesheet read in from --input
    ch_dataset            // channel: datasets with sample names and paths
    ch_batches            // channel: batch file information
    ch_classes            // channel: classes from samplesheet metadata
    ch_annotations        // channel: gene annotations from the GEX_REF_PREPROCESSOR subworkflow
    gex_genes_file        // classifier-selected GEX feature list, or null
    gex_norm_factors_file  // classifier-trained normalization factors, or null
    class_colors_file      // classifier class colors, or null
    random_seed           //  int: Random seed for reproducibility

    main:

    ch_versions = channel.empty()
    ch_visuals = channel.empty()
    ch_norm_factors = channel.empty()
    def gex_matrix_for_normalization = null

    CONCATENATE_GEX (
        ch_dataset
    )
    ch_versions = ch_versions.mix(CONCATENATE_GEX.out.versions)

    FILTER_BY_ANNOTATION (
        CONCATENATE_GEX.out.dataset_name,
        CONCATENATE_GEX.out.concatenated_gex_csv,
        ch_annotations,
    )
    ch_versions = ch_versions.mix(FILTER_BY_ANNOTATION.out.versions)

    MERGE_DATASETS (
        FILTER_BY_ANNOTATION.out.filtered_annotation_csv.collect()
            .map { filtered_genes_csvs -> filtered_genes_csvs.sort { a, b -> a.name <=> b.name } }
    )
    ch_versions = ch_versions.mix(MERGE_DATASETS.out.versions)
    gex_matrix_for_normalization = MERGE_DATASETS.out.merged_csv

    if (params.gex_intermediate_umaps && !params.skip_umaps) {
        if (params.tsne) {
            TSNE_GEX_RAW_BY_BATCH (
                "gex.raw.by_batch",
                MERGE_DATASETS.out.merged_csv,
                ch_batches,
                false,
                class_colors_file ?: '',
                random_seed
            )
            ch_versions = ch_versions.mix(TSNE_GEX_RAW_BY_BATCH.out.versions)
            ch_visuals = ch_visuals.mix(TSNE_GEX_RAW_BY_BATCH.out.tsne_svg)

            TSNE_GEX_RAW_BY_CLASS (
                "gex.raw.by_class",
                MERGE_DATASETS.out.merged_csv,
                ch_classes,
                false,
                class_colors_file ?: '',
                random_seed
            )
            ch_versions = ch_versions.mix(TSNE_GEX_RAW_BY_CLASS.out.versions)
            ch_visuals = ch_visuals.mix(TSNE_GEX_RAW_BY_CLASS.out.tsne_svg)
        }
        else {
            UMAP_GEX_RAW_BY_BATCH (
                "gex.raw.by_batch",
                MERGE_DATASETS.out.merged_csv,
                ch_batches,
                false,
                class_colors_file ?: '',
                random_seed
            )
            ch_versions = ch_versions.mix(UMAP_GEX_RAW_BY_BATCH.out.versions)
            ch_visuals = ch_visuals.mix(UMAP_GEX_RAW_BY_BATCH.out.umap_svg)

            UMAP_GEX_RAW_BY_CLASS (
                "gex.raw.by_class",
                MERGE_DATASETS.out.merged_csv,
                ch_classes,
                false,
                class_colors_file ?: '',
                random_seed
            )
            ch_versions = ch_versions.mix(UMAP_GEX_RAW_BY_CLASS.out.versions)
            ch_visuals = ch_visuals.mix(UMAP_GEX_RAW_BY_CLASS.out.umap_svg)
        }
    }

    if (params.batch_correction) {
        BATCH_CORRECT (
            MERGE_DATASETS.out.merged_csv,
            ch_batches,
            ch_classes
        )
        ch_versions = ch_versions.mix(BATCH_CORRECT.out.versions)
        gex_matrix_for_normalization = BATCH_CORRECT.out.batch_corrected_csv

        if (params.gex_intermediate_umaps && !params.skip_umaps) {
            if (params.tsne) {
                TSNE_GEX_BC_BY_BATCH (
                    "gex.batch_corrected.by_batch",
                    BATCH_CORRECT.out.batch_corrected_csv,
                    ch_batches,
                    false,
                    class_colors_file ?: '',
                    random_seed
                )
                ch_versions = ch_versions.mix(TSNE_GEX_BC_BY_BATCH.out.versions)
                ch_visuals = ch_visuals.mix(TSNE_GEX_BC_BY_BATCH.out.tsne_svg)

                TSNE_GEX_BC_BY_CLASS (
                    "gex.batch_corrected.by_class",
                    BATCH_CORRECT.out.batch_corrected_csv,
                    ch_classes,
                    false,
                    class_colors_file ?: '',
                    random_seed
                )
                ch_versions = ch_versions.mix(TSNE_GEX_BC_BY_CLASS.out.versions)
                ch_visuals = ch_visuals.mix(TSNE_GEX_BC_BY_CLASS.out.tsne_svg)
            }
            else {
                UMAP_GEX_BC_BY_BATCH (
                    "gex.batch_corrected.by_batch",
                    BATCH_CORRECT.out.batch_corrected_csv,
                    ch_batches,
                    false,
                    class_colors_file ?: '',
                    random_seed
                )
                ch_versions = ch_versions.mix(UMAP_GEX_BC_BY_BATCH.out.versions)
                ch_visuals = ch_visuals.mix(UMAP_GEX_BC_BY_BATCH.out.umap_svg)

                UMAP_GEX_BC_BY_CLASS (
                    "gex.batch_corrected.by_class",
                    BATCH_CORRECT.out.batch_corrected_csv,
                    ch_classes,
                    false,
                    class_colors_file ?: '',
                    random_seed
                )
                ch_versions = ch_versions.mix(UMAP_GEX_BC_BY_CLASS.out.versions)
                ch_visuals = ch_visuals.mix(UMAP_GEX_BC_BY_CLASS.out.umap_svg)
            }
        }
    }

    NORMALIZE (
        gex_matrix_for_normalization,
        ch_annotations,
        gex_norm_factors_file
    )
    ch_versions = ch_versions.mix(NORMALIZE.out.versions)
    ch_norm_factors = NORMALIZE.out.norm_factors_rds

    def gex_processed_matrix

    if (gex_genes_file) {
        FILTER_BY_LIST (
            "merged_datasets",
            NORMALIZE.out.normalized_csv,
            gex_genes_file
        )
        ch_versions = ch_versions.mix(FILTER_BY_LIST.out.versions)
        gex_processed_matrix = FILTER_BY_LIST.out.filtered_by_list_csv.map { f -> [ 'merged_datasets', 'merged_datasets', f ] }
    }
    else {
        gex_processed_matrix = NORMALIZE.out.normalized_csv.map { f -> [ 'merged_datasets', 'merged_datasets', f ] }
    }

    if (params.gex_intermediate_umaps && !params.skip_umaps) {
        if (params.tsne) {
            TSNE_GEX_NORM_BY_BATCH (
                "gex.normalized.by_batch",
                NORMALIZE.out.normalized_csv,
                ch_batches,
                false,
                class_colors_file ?: '',
                random_seed
            )
            ch_versions = ch_versions.mix(TSNE_GEX_NORM_BY_BATCH.out.versions)
            ch_visuals = ch_visuals.mix(TSNE_GEX_NORM_BY_BATCH.out.tsne_svg)

            TSNE_GEX_NORM_BY_CLASS (
                "gex.normalized.by_class",
                NORMALIZE.out.normalized_csv,
                ch_classes,
                false,
                class_colors_file ?: '',
                random_seed
            )
            ch_versions = ch_versions.mix(TSNE_GEX_NORM_BY_CLASS.out.versions)
            ch_visuals = ch_visuals.mix(TSNE_GEX_NORM_BY_CLASS.out.tsne_svg)
        }
        else {
            UMAP_GEX_NORM_BY_BATCH (
                "gex.normalized.by_batch",
                NORMALIZE.out.normalized_csv,
                ch_batches,
                false,
                class_colors_file ?: '',
                random_seed
            )
            ch_versions = ch_versions.mix(UMAP_GEX_NORM_BY_BATCH.out.versions)
            ch_visuals = ch_visuals.mix(UMAP_GEX_NORM_BY_BATCH.out.umap_svg)

            UMAP_GEX_NORM_BY_CLASS (
                "gex.normalized.by_class",
                NORMALIZE.out.normalized_csv,
                ch_classes,
                false,
                class_colors_file ?: '',
                random_seed
            )
            ch_versions = ch_versions.mix(UMAP_GEX_NORM_BY_CLASS.out.versions)
            ch_visuals = ch_visuals.mix(UMAP_GEX_NORM_BY_CLASS.out.umap_svg)
        }
    }

    if (params.tsne) {
        TSNE_GEX_PROCESSED_BY_BATCH (
            "gex.processed.by_batch",
            gex_processed_matrix.map { _dataset_name, _sample_name, f -> f },
            ch_batches,
            false,
            class_colors_file ?: '',
            random_seed
        )
        ch_versions = ch_versions.mix(TSNE_GEX_PROCESSED_BY_BATCH.out.versions)
        ch_visuals = ch_visuals.mix(TSNE_GEX_PROCESSED_BY_BATCH.out.tsne_svg)

        TSNE_GEX_PROCESSED_BY_CLASS (
            "gex.processed.by_class",
            gex_processed_matrix.map { _dataset_name, _sample_name, f -> f },
            ch_classes,
            false,
            class_colors_file ?: '',
            random_seed
        )
        ch_versions = ch_versions.mix(TSNE_GEX_PROCESSED_BY_CLASS.out.versions)
        ch_visuals = ch_visuals.mix(TSNE_GEX_PROCESSED_BY_CLASS.out.tsne_svg)
    }
    else if (!params.skip_umaps) {
        UMAP_GEX_PROCESSED_BY_BATCH (
            "gex.processed.by_batch",
            gex_processed_matrix.map { _dataset_name, _sample_name, f -> f },
            ch_batches,
            false,
            class_colors_file ?: '',
            random_seed
        )
        ch_versions = ch_versions.mix(UMAP_GEX_PROCESSED_BY_BATCH.out.versions)
        ch_visuals = ch_visuals.mix(UMAP_GEX_PROCESSED_BY_BATCH.out.umap_svg)

        UMAP_GEX_PROCESSED_BY_CLASS (
            "gex.processed.by_class",
            gex_processed_matrix.map { _dataset_name, _sample_name, f -> f },
            ch_classes,
            false,
            class_colors_file ?: '',
            random_seed
        )
        ch_versions = ch_versions.mix(UMAP_GEX_PROCESSED_BY_CLASS.out.versions)
        ch_visuals = ch_visuals.mix(UMAP_GEX_PROCESSED_BY_CLASS.out.umap_svg)
    }

    TRANSPOSE (
        gex_processed_matrix.map { _dataset_name, _sample_name, f -> f },
        "gex"
    )
    ch_versions = ch_versions.mix(TRANSPOSE.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'gexml_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        )


    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    transposed_csv = TRANSPOSE.out.transposed_csv
    norm_factors_rds = ch_norm_factors
    visuals        = ch_visuals                  // channel: path(*.png|*.svg)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
