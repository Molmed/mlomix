/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPROCESS_MINFI       } from '../modules/local/dnam/preprocess_minfi/main'
include { CONCATENATE_DNAM       } from '../modules/local/dnam/concatenate_dnam'
include { MERGE_DATASETS         } from '../modules/local/dnam/merge_datasets'
include { COMPRESS_DNAM          } from '../modules/local/dnam/compress_dnam/main'
include { P_VAL_CORRECTION       } from '../modules/local/dnam/p_val_correction/main'
include { FILTER_PROBES_BY_LIST  } from '../modules/local/dnam/filter_probes_by_list/main'
include { FILTER_BY_MISSING      } from '../modules/local/dnam/filter_by_missing/main'
include { FILTER_BY_VARIANCE     } from '../modules/local/filter_by_variance/main'
include { IMPUTE                 } from '../modules/local/dnam/impute/main'
include { UMAP as UMAP_DNAM_BY_CLASS } from '../modules/local/umap/umap'
include { TSNE as TSNE_DNAM_BY_CLASS } from '../modules/local/tsne/tsne'
include { TRANSPOSE              } from '../modules/local/transpose'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DNAM {

    take:
    ch_samplesheet
    ch_classes
    random_seed
    _ch_precomputed_beta_matrix
    _ch_precomputed_pvals
    _ch_use_precomputed_dnam

    main:

    ch_versions = channel.empty()
    ch_visuals = channel.empty()
    ch_dnam_probe_list = channel.value(file(params.dnam_probes_file ?: params.common_probes, checkIfExists: true))
    def classifier_configured = (params.classifier_name != null) && (params.classifier_version != null)
    def skip_missing_filter = (params.dnam_probes_file != null) || params.dnam_skip_filter_by_missing || classifier_configured
    def skip_variance_filter = (params.dnam_probes_file != null) || params.dnam_skip_filter_by_variance || classifier_configured

    //
    // Build per-sample DNAm channels from precomputed inputs
    //
    ch_precomputed_rows = ch_samplesheet
        .flatMap { rows -> (rows instanceof List) ? rows : [rows] }
        .filter { row -> row['dnam_beta_matrix_file'] }

    ch_precomputed_pairs = ch_precomputed_rows
        .filter { row -> row['dnam_pvals_file'] }
        .map { row ->
            def dataset_name = row['dataset'] == null ? 'dataset' : row['dataset'].toString().trim()
            dataset_name = dataset_name ? dataset_name : 'dataset'
            dataset_name = dataset_name.replaceAll(/[^A-Za-z0-9._-]/, '_')
            def sample_name = row['id'] == null ? 'sample' : row['id'].toString().trim()
            sample_name = sample_name ? sample_name : 'sample'
            sample_name = sample_name.replaceAll(/[^A-Za-z0-9._-]/, '_')
            [
                dataset_name,
                sample_name,
                file(row['dnam_beta_matrix_file'], checkIfExists: true),
                file(row['dnam_pvals_file'], checkIfExists: true)
            ]
        }

    ch_precomputed_beta_only = ch_precomputed_rows
        .filter { row -> !row['dnam_pvals_file'] }
        .map { row ->
            def dataset_name = row['dataset'] == null ? 'dataset' : row['dataset'].toString().trim()
            dataset_name = dataset_name ? dataset_name : 'dataset'
            dataset_name = dataset_name.replaceAll(/[^A-Za-z0-9._-]/, '_')
            def sample_name = row['id'] == null ? 'sample' : row['id'].toString().trim()
            sample_name = sample_name ? sample_name : 'sample'
            sample_name = sample_name.replaceAll(/[^A-Za-z0-9._-]/, '_')
            [
                dataset_name,
                sample_name,
                file(row['dnam_beta_matrix_file'], checkIfExists: true)
            ]
        }

    //
    // MODULE: Preprocess methylation array data with minfi
    //
    ch_minfi_samples = ch_samplesheet
        .flatMap { rows -> (rows instanceof List) ? rows : [rows] }
        .filter { row -> row['sentrix_id'] || row['sentrix_position'] || row['idats_basename'] }
        .map { row ->
            def dataset_name = row['dataset'] == null ? 'dataset' : row['dataset'].toString().trim()
            dataset_name = dataset_name ? dataset_name : 'dataset'
            dataset_name = dataset_name.replaceAll(/[^A-Za-z0-9._-]/, '_')
            def sample_name = row['id'] == null ? 'sample' : row['id'].toString().trim()
            sample_name = sample_name ? sample_name : 'sample'
            sample_name = sample_name.replaceAll(/[^A-Za-z0-9._-]/, '_')
            [
                dataset_name,
                sample_name,
                row['sentrix_id'],
                row['sentrix_position'],
                row['idats_basename']
            ]
        }

    PREPROCESS_MINFI (
        ch_minfi_samples
    )
    ch_versions = ch_versions.mix(PREPROCESS_MINFI.out.versions)

    ch_minfi_pairs = PREPROCESS_MINFI.out.corrected_inputs

    ch_pairs_with_pvals = ch_minfi_pairs.mix(ch_precomputed_pairs)

    ch_betas_for_filtering = ch_pairs_with_pvals
        .map { dataset_name, sample_name, betas, _detection_pvals -> [dataset_name, sample_name, betas] }
        .mix(ch_precomputed_beta_only)

    //
    // MODULE: Filter DNAm betas per sample (common probes)
    //
    FILTER_PROBES_BY_LIST (
        ch_betas_for_filtering,
        ch_dnam_probe_list
    )
    ch_versions = ch_versions.mix(FILTER_PROBES_BY_LIST.out.versions)

    ch_filtered_betas_for_compression = FILTER_PROBES_BY_LIST.out.filtered_betas

    COMPRESS_DNAM (
        ch_filtered_betas_for_compression
    )
    ch_versions = ch_versions.mix(COMPRESS_DNAM.out.versions)

    ch_compressed_by_key = COMPRESS_DNAM.out.compressed_betas
        .map { dataset_name, sample_name, compressed_betas ->
            ["${dataset_name}__${sample_name}", dataset_name, sample_name, compressed_betas]
        }

    ch_pvals_by_key = ch_pairs_with_pvals
        .map { dataset_name, sample_name, _betas, detection_pvals ->
            ["${dataset_name}__${sample_name}", detection_pvals]
        }

    ch_compressed_with_optional_pvals = ch_compressed_by_key
        .join(ch_pvals_by_key, remainder: true)

    ch_pairs_for_correction = ch_compressed_with_optional_pvals
        .filter { tuple_item -> tuple_item[4] != null }
        .map { tuple_item -> [tuple_item[1], tuple_item[2], tuple_item[3], tuple_item[4]] }

    ch_compressed_beta_only = ch_compressed_with_optional_pvals
        .filter { tuple_item -> tuple_item[4] == null }
        .map { tuple_item -> [tuple_item[1], tuple_item[2], tuple_item[3]] }

    //
    // MODULE: Replace beta values with NaN where detection p-value >= threshold
    // Run per sample, then concatenate corrected betas by dataset.
    //
    P_VAL_CORRECTION (
        ch_pairs_for_correction
    )
    ch_versions = ch_versions.mix(P_VAL_CORRECTION.out.versions)

    ch_corrected_or_passthrough_betas = P_VAL_CORRECTION.out.corrected_betas
        .mix(ch_compressed_beta_only)

    ch_common_filtered_by_dataset = ch_corrected_or_passthrough_betas
        .groupTuple()
        .map { dataset_name, sample_names, beta_paths ->
            // Sort by sample_name to ensure deterministic chunk order for caching
            def sorted_pairs = [sample_names, beta_paths].transpose().sort { a, b -> a[0] <=> b[0] }
            def sorted_samples = sorted_pairs.collect { pair -> pair[0] }
            def sorted_paths = sorted_pairs.collect { pair -> pair[1] }
            [
                dataset_name,
                sorted_samples,
                sorted_paths
            ]
        }
        .flatMap { dataset_name, sample_names, beta_paths ->
            def chunk_size = Math.max(1, (params.dnam_samples_per_chunk ?: 200) as Integer)
            def sample_chunks = sample_names.collate(chunk_size)
            def path_chunks = beta_paths.collate(chunk_size)

            sample_chunks.indices.collect { idx ->
                [dataset_name, idx, sample_chunks[idx], path_chunks[idx]]
            }
        }

    CONCATENATE_DNAM (
        ch_common_filtered_by_dataset
    )
    ch_versions = ch_versions.mix(CONCATENATE_DNAM.out.versions)

    MERGE_DATASETS (
        CONCATENATE_DNAM.out.beta_matrix.collect()
            .map { beta_matrices -> beta_matrices.sort { a, b -> a.name <=> b.name } }
    )
    ch_versions = ch_versions.mix(MERGE_DATASETS.out.versions)

    ch_dataset_betas_for_missing = MERGE_DATASETS.out.beta_matrix
        .map { beta_matrix ->
            ["merged_datasets", "merged_datasets", beta_matrix]
        }

    ch_missing_filtered_betas = ch_dataset_betas_for_missing
    if (!skip_missing_filter) {
        FILTER_BY_MISSING (
            ch_dataset_betas_for_missing
        )
        ch_versions = ch_versions.mix(FILTER_BY_MISSING.out.versions)
        ch_missing_filtered_betas = FILTER_BY_MISSING.out.missing_filtered_betas
    }

    ch_variance_input_betas = ch_missing_filtered_betas
    if (!skip_variance_filter) {
        FILTER_BY_VARIANCE (
            ch_missing_filtered_betas,
            'dnam'
        )
        ch_versions = ch_versions.mix(FILTER_BY_VARIANCE.out.versions)
        ch_variance_input_betas = FILTER_BY_VARIANCE.out.variance_filtered_betas

        ch_visuals = ch_visuals.mix(FILTER_BY_VARIANCE.out.variance_plot_before_svg)
        ch_visuals = ch_visuals.mix(FILTER_BY_VARIANCE.out.variance_plot_after_svg)
    }

    IMPUTE (
        ch_variance_input_betas
    )
    ch_versions = ch_versions.mix(IMPUTE.out.versions)

    ch_imputed_betas = IMPUTE.out.imputed_betas

    if (params.tsne) {
        TSNE_DNAM_BY_CLASS (
            "dnam.by_class",
            ch_imputed_betas.map { _dataset_name, _sample_name, betas -> betas },
            ch_classes,
            false,
            params.class_colors_file ?: '',
            random_seed
        )
        ch_versions = ch_versions.mix(TSNE_DNAM_BY_CLASS.out.versions)
        ch_visuals = ch_visuals.mix(TSNE_DNAM_BY_CLASS.out.tsne_svg)
    }
    else if (!params.skip_umaps) {
        UMAP_DNAM_BY_CLASS (
            "dnam.by_class",
            ch_imputed_betas.map { _dataset_name, _sample_name, betas -> betas },
            ch_classes,
            false,
            params.class_colors_file ?: '',
            random_seed
        )
        ch_versions = ch_versions.mix(UMAP_DNAM_BY_CLASS.out.versions)
        ch_visuals = ch_visuals.mix(UMAP_DNAM_BY_CLASS.out.umap_svg)
    }

    TRANSPOSE (
        ch_imputed_betas.map { _dataset_name, _sample_name, betas -> betas },
        "dnam"
    )
    ch_versions = ch_versions.mix(TRANSPOSE.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'methylml_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    versions = ch_versions // channel: [ path(versions.yml) ]
    transposed_csv = TRANSPOSE.out.transposed_csv
    visuals = ch_visuals // channel: path(*.png|*.svg)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
