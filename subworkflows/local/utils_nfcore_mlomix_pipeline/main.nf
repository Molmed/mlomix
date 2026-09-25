//
// Subworkflow with functionality specific to the nf-core/mlomix pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { CLASS_FILTER_AND_REPORT   } from '../../../modules/local/mlomix/class_reporter/main'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

def resolveClassifierConfigFile(manifest, configDir, key) {
    def relativePath = manifest[key]?.toString()
    if (!relativePath) {
        return null
    }

    def resolvedPath = new File(configDir, relativePath).canonicalFile
    if (!resolvedPath.path.startsWith(configDir.path + File.separator)) {
        error("Classifier config file escapes classifier directory: ${key}")
    }
    if (!resolvedPath.isFile()) {
        error("Classifier config file not found for '${key}': ${resolvedPath}")
    }
    return resolvedPath.path
}

def resolveClassifierConfig() {
    if (!params.classifier_name && !params.classifier_version) {
        return [
            genome: params.genome,
            annotation_version: params.annotation_version ?: 109,
            gex_genes_file: params.gex_genes_file,
            gex_norm_factors_file: params.gex_norm_factors_file,
            dnam_probes_file: params.dnam_probes_file,
            dnam_classes_file: params.dnam_classes_file,
            class_colors_file: params.class_colors_file,
        ]
    }

    if (!params.classifier_name || !params.classifier_version) {
        error("Both --classifier-name and --classifier-version must be provided together.")
    }

    def classifierName = params.classifier_name.toString()
    def classifierVersion = params.classifier_version.toString()
    def safePathPart = /^[A-Za-z0-9_.-]+$/
    if (!(classifierName ==~ safePathPart) || !(classifierVersion ==~ safePathPart)) {
        error("Classifier name and version may contain only letters, numbers, '.', '_' and '-'.")
    }

    def classifierRoot = file(params.classifier_configs_dir).toFile().canonicalFile
    def configDir = new File(classifierRoot, "${classifierName}/${classifierVersion}").canonicalFile
    if (!configDir.isDirectory()) {
        configDir = new File(classifierRoot, "${classifierName}/v${classifierVersion}").canonicalFile
    }
    if (!configDir.isDirectory()) {
        error("Classifier configuration directory not found for ${classifierName} ${classifierVersion} under ${classifierRoot}")
    }

    def manifestPath = new File(configDir, 'manifest.yml')
    if (!manifestPath.isFile()) {
        error("Classifier manifest not found: ${manifestPath}")
    }

    def manifest = new org.yaml.snakeyaml.Yaml().load(manifestPath.text)
    if (!(manifest instanceof Map) || !(manifest.classifier instanceof Map)) {
        error("Classifier manifest is missing the 'classifier' section: ${manifestPath}")
    }
    if (manifest.classifier.name?.toString() != classifierName ||
        manifest.classifier.version?.toString() != classifierVersion) {
        error("Classifier manifest identity does not match its directory: ${manifestPath}")
    }

    def resolved = [
        genome: params.genome ?: manifest['genome'],
        annotation_version: params.annotation_version ?: manifest['annotation_version'] ?: 109,
        gex_genes_file: params.gex_genes_file,
        gex_norm_factors_file: params.gex_norm_factors_file,
        dnam_probes_file: params.dnam_probes_file,
        dnam_classes_file: params.dnam_classes_file,
        class_colors_file: params.class_colors_file,
    ]
    if (!resolved.gex_genes_file) {
        resolved.gex_genes_file = resolveClassifierConfigFile(manifest, configDir, 'gex_genes_file')
    }
    if (!resolved.gex_norm_factors_file) {
        resolved.gex_norm_factors_file = resolveClassifierConfigFile(manifest, configDir, 'gex_norm_factors_file')
    }
    if (!resolved.dnam_probes_file) {
        resolved.dnam_probes_file = resolveClassifierConfigFile(manifest, configDir, 'dnam_probes_file')
    }
    if (!resolved.dnam_classes_file) {
        resolved.dnam_classes_file = resolveClassifierConfigFile(manifest, configDir, 'dnam_classes_file')
    }
    if (!resolved.class_colors_file) {
        resolved.class_colors_file = resolveClassifierConfigFile(manifest, configDir, 'class_colors_file')
    }
    log.info("Using classifier configuration ${classifierName} ${classifierVersion} (annotation_version=${resolved.annotation_version})")
    return resolved
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs  // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    _input            //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message
    classifier_config // map: Resolved classifier configuration

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    before_text = """
-\033[2m----------------------------------------------------\033[0m-
                                        \033[0;32m,--.\033[0;30m/\033[0;32m,-.\033[0m
\033[0;34m        ___     __   __   __   ___     \033[0;32m/,-._.--~\'\033[0m
\033[0;34m  |\\ | |__  __ /  ` /  \\ |__) |__         \033[0;33m}  {\033[0m
\033[0;34m  | \\| |       \\__, \\__/ |  \\ |___     \033[0;32m\\`-._,-`-,\033[0m
                                        \033[0;32m`._,._,\'\033[0m
\033[0;35m  nf-core/mlomix ${workflow.manifest.version}\033[0m
-\033[2m----------------------------------------------------\033[0m-
"""
    after_text = """${workflow.manifest.doi ? "\n* The pipeline\n" : ""}${workflow.manifest.doi.tokenize(",").collect { doi -> "    https://doi.org/${doi.trim().replace('https://doi.org/','')}"}.join("\n")}${workflow.manifest.doi ? "\n" : ""}
* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/mlomix/blob/main/CITATIONS.md
"""
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Shared class report (filter + report), generated once before any workflow branching
    //
    ch_samplesheet_report = channel.fromPath(params.input, checkIfExists: true)
    CLASS_FILTER_AND_REPORT(
        ch_samplesheet_report,
        classifier_config.gex_classes_file ?: '',
        classifier_config.dnam_classes_file ?: ''
    )
    ch_versions = ch_versions.mix(CLASS_FILTER_AND_REPORT.out.versions)
    ch_class_plot_gex_filtered = CLASS_FILTER_AND_REPORT.out.class_plot_gex_filtered
    ch_class_plot_dnam_filtered = CLASS_FILTER_AND_REPORT.out.class_plot_dnam_filtered
    ch_classes_gex_filtered  = CLASS_FILTER_AND_REPORT.out.classes_gex_filtered
    ch_classes_dnam_filtered = CLASS_FILTER_AND_REPORT.out.classes_dnam_filtered
    ch_filtered_samplesheet_csv = CLASS_FILTER_AND_REPORT.out.samplesheet_filtered_csv

    ch_rows = ch_filtered_samplesheet_csv
        .flatMap { csv_file ->
            samplesheetToList(csv_file, "${projectDir}/assets/schema_input.json")
        }
        .map { row -> row[0] }
        .collect()

    ch_mode_info = ch_rows
        .map { filtered_rows ->
            validateInputSamplesheetModes(filtered_rows)
        }
        .map { mode_info ->
            log.info "Parsed samplesheet modes: run_gex=${mode_info.run_gex}, run_dnam=${mode_info.run_dnam}, use_precomputed_dnam=${mode_info.use_precomputed_dnam}"

            if (mode_info.run_gex && !classifier_config.genome) {
                error("GEX input was detected but '--genome' is missing. Please provide a supported genome key (e.g. --genome GRCh38).")
            }

            if (mode_info.run_gex && !classifier_config.annotation_version) {
                error("GEX input was detected but '--annotation_version' is missing. Please provide an Ensembl release (e.g. --annotation_version 109).")
            }

            return mode_info
        }

    ch_run_gex = ch_mode_info.map { it.run_gex }
    ch_run_dnam = ch_mode_info.map { it.run_dnam }
    ch_use_precomputed_dnam = ch_mode_info.map { it.use_precomputed_dnam }

    ch_rows
        .flatMap { rows -> rows }
        .filter { row -> row.gex_feature_counts_file }
        .set { ch_gex_rows }

    // Group GEX rows by dataset to build CONCATENATE_GEX input tuples.
    ch_gex_rows
        .map { row -> [ row.dataset, row ] }
        .groupTuple()
        .set { ch_gex_samplesheet }

    ch_gex_samplesheet
        .map { dataset ->
            def sample_names = dataset[1].collect { sample -> sample.id }
            def sample_paths = dataset[1].collect { sample -> file(sample.gex_feature_counts_file, checkIfExists: true) }
            [dataset[0], sample_names, sample_paths]
        }
        .set { ch_datasets }

    ch_gex_samplesheet
        .map { dataset ->
            dataset[1].collect { sample ->
                def batch_name = sample.batch ?: ''
                def suffix = batch_name ? "_${batch_name}" : ''
                "${sample.id}\t${dataset[0]}${suffix}"
            }.join('\n')
        }
        .map { content -> "sample\tbatch\n${content}\n" }
        .collectFile(name: "batches.tsv",
                     newLine: false,
                     keepHeader: true,
                     storeDir: "${params.outdir}/batch")
        .set { ch_batches }

    ch_rows
        .map { grouped_rows ->
            def rows = grouped_rows ?: []
            def header = "sample\tclass"
            def body = rows.collect { data ->
                def sample_name = data['sample'] ?: data['id']
                "${sample_name}\t${data['class'] ?: ''}"
            }.join('\n')
            "${header}\n${body}\n"
        }
        .collectFile(name: 'classes.tsv',
                     newLine: false,
                     storeDir: "${params.outdir}/class")
        .set { ch_classes }

    // Create separate channel for DNAM using the full rows (same data as GEX, pre-filtered by CLASS_FILTER_AND_REPORT separately)
    ch_rows
        .set { ch_dnam_samplesheet }

    ch_dnam_beta_matrix = channel.empty()
    ch_dnam_pvals = channel.empty()
    ch_genome = channel.value(classifier_config.genome)
    ch_annotation_version = channel.value(classifier_config.annotation_version)
    ch_random_seed = channel.value(params.random_seed)

    emit:
    gex_samplesheet       = ch_gex_samplesheet
    datasets              = ch_datasets
    batches               = ch_batches
    classes               = ch_classes
    dnam_samplesheet      = ch_dnam_samplesheet
    dnam_beta_matrix      = ch_dnam_beta_matrix
    dnam_pvals            = ch_dnam_pvals
    annotation_version    = ch_annotation_version
    genome                = ch_genome
    random_seed           = ch_random_seed
    run_gex               = ch_run_gex
    run_dnam              = ch_run_dnam
    use_precomputed_dnam  = ch_use_precomputed_dnam
    classes_gex_filtered  = ch_classes_gex_filtered
    classes_dnam_filtered = ch_classes_dnam_filtered
    class_plot_gex_filtered = ch_class_plot_gex_filtered
    class_plot_dnam_filtered = ch_class_plot_dnam_filtered
    samplesheet_filtered_csv = ch_filtered_samplesheet_csv
    versions              = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate modes and optional columns in combined input samplesheet
//
def validateInputSamplesheetModes(rows) {
    if (!rows || rows.isEmpty()) {
        error('Input samplesheet is empty after parsing.')
    }

    def run_gex = rows.any { row -> row.gex_feature_counts_file }

    // A sample can carry both modalities; DNAM detection must be based on DNAM columns, not on absence of GEX.
    def dnam_rows = rows.findAll { row ->
        row['dnam_beta_matrix_file'] || row['dnam_pvals_file'] || row['sentrix_id'] || row['sentrix_position'] || row['idats_basename']
    }
    def run_dnam = !dnam_rows.isEmpty()

    // Validate that each DNAM sample uses exactly one mode (per-sample validation)
    dnam_rows.each { row ->
        def has_precomputed = row['dnam_beta_matrix_file'] || row['dnam_pvals_file']
        def has_idat = row['sentrix_id'] || row['sentrix_position'] || row['idats_basename']

        // Check for incomplete precomputed mode
        if (row['dnam_pvals_file'] && !row['dnam_beta_matrix_file']) {
            error("DNAM sample '${row.id}' has invalid precomputed mode: dnam_pvals_file requires dnam_beta_matrix_file.")
        }
        if (!row['dnam_beta_matrix_file'] && !has_idat) {
            error("DNAM sample '${row.id}' has no valid DNAM input. Provide either dnam_beta_matrix_file (optionally with dnam_pvals_file) or (sentrix_id + sentrix_position + idats_basename).")
        }

        // Check for incomplete IDAT mode
        if ((row['sentrix_id'] || row['sentrix_position'] || row['idats_basename']) &&
            (!row['sentrix_id'] || !row['sentrix_position'] || !row['idats_basename'])) {
            error("DNAM sample '${row.id}' has incomplete IDAT mode: sentrix_id, sentrix_position and idats_basename are all required together.")
        }

        // Check for mixing modes within a single sample
        if (has_precomputed && has_idat) {
            error("DNAM sample '${row.id}' mixes precomputed and IDAT modes. Each sample must use exactly one mode.")
        }

        // Check that sample has at least one DNAM mode
        if (!has_precomputed && !has_idat) {
            error("DNAM sample '${row.id}' has no valid DNAM input. Provide either dnam_beta_matrix_file (optionally with dnam_pvals_file) or (sentrix_id + sentrix_position + idats_basename).")
        }
    }

    // Mixed-mode support: allow different samples to use different DNAM modes in the same run
    def precomputed_rows = dnam_rows.findAll { row -> row['dnam_beta_matrix_file'] || row['dnam_pvals_file'] }

    def use_precomputed_dnam = !precomputed_rows.isEmpty()

    if (!run_gex && !run_dnam) {
        error('No runnable rows found. Provide gex_feature_counts_file and/or valid DNAM input columns.')
    }

    [
        run_gex: run_gex,
        run_dnam: run_dnam,
        use_precomputed_dnam: use_precomputed_dnam,
        precomputed_rows: precomputed_rows,
        idat_rows: dnam_rows.findAll { row -> row['sentrix_id'] || row['sentrix_position'] || row['idats_basename'] }
    ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
