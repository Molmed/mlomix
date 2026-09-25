process NORMALIZE {
    label 'process_medium'

    conda "bioconda::bioconductor-edger=4.4.0"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bioconductor-edger:4.4.0--r44h3df3fcb_0 '
        : 'biocontainers/bioconductor-edger:4.4.0--r44h3df3fcb_0 '}"

    input:
    path merged_data_path
    path ref_path
    val gex_norm_factors_file

    output:
    path "normalized.csv", emit: normalized_csv
    path "gex_norm_factors.rds", emit: norm_factors_rds
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def norm_factors_file = gex_norm_factors_file ?: ''

    """
    #!/usr/bin/env Rscript

    library(edgeR)

    normalize <- function(gex_path, ref_path, output_path) {
        x <- read.csv(gex_path, row.names = 1, header= TRUE, check.names = FALSE)
        annot <- read.csv(ref_path, row.names = 1, header= TRUE, check.names = FALSE)

        x_length_norm <- (x * 10^3) / annot\$length
        factors_path <- "gex_norm_factors.rds"

        d <- DGEList(counts = x_length_norm)

        if (nzchar("${norm_factors_file}")) {
            message(sprintf("NORMALIZE: Using provided gex norm factors file: %s", "${norm_factors_file}"))
            train_params <- readRDS("${norm_factors_file}")

            # Use the training norm factor median for all test samples.
            d\$samples\$norm.factors <- rep(median(train_params\$norm_factors), ncol(x))

            CPM <- cpm(d, log = TRUE)
            write.csv(CPM, output_path)
            saveRDS(train_params, factors_path)
        } else {
            message("NORMALIZE: No gex_norm_factors_file provided; computing TMM normalization factors from current input")
            TMM <- calcNormFactors(d, method = "TMM")
            saveRDS(list(
                norm_factors = TMM\$samples\$norm.factors,
                ref_lib_size = mean(TMM\$samples\$lib.size)
            ), factors_path)

            CPM <- cpm(TMM, log = TRUE)
            write.csv(CPM, output_path)
        }
    }

    normalize("${merged_data_path}", "${ref_path}", "normalized.csv")

    # Create versions file
    writeLines(c(
        '"${task.process}":',
        paste0('    R: "', R.version.string, '"'),
        paste0('    E: "', packageVersion("edgeR"), '"')
    ), "versions.yml")
    """
}
