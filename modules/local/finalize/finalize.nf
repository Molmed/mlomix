process FINALIZE {
    label 'process_high'

    conda "conda-forge::pandas=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.3' :
        'biocontainers/pandas:2.2.3' }"

    input:
    path classes_file
    path dnam_transposed
    path gex_transposed
    path gex_norm_factors_rds
    path visual_files
    path class_plot_gex_filtered
    path class_plot_dnam_filtered
    path gex_classes_filtered
    path dnam_classes_filtered
    path samplesheet_filtered_csv
    path prediction_files

    output:
    path "labels.csv", emit: labels_csv
    path "features.dnam.csv", emit: dnam_features_csv
    path "features.gex.csv", emit: gex_features_csv
    path "feature_names.dnam.txt", emit: dnam_feature_names
    path "feature_names.gex.txt", emit: gex_feature_names
    path "gex_norm_factors.rds", emit: gex_norm_factors_rds
    path "class_distribution.gex.filtered.svg", emit: class_distribution_gex_filtered_svg
    path "class_distribution.dnam.filtered.svg", emit: class_distribution_dnam_filtered_svg
    path "classes.gex.filtered.tsv", emit: classes_gex_filtered
    path "classes.dnam.filtered.tsv", emit: classes_dnam_filtered
    path "samplesheet.filtered.csv", emit: samplesheet_filtered_csv
    path "visuals", emit: visuals_dir
    path "predictions.*.csv", optional: true, emit: predictions
    path "class_names.gex.txt", emit: class_names_gex
    path "class_names.dnam.txt", emit: class_names_dnam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def visual_file_names = (visual_files instanceof List ? visual_files : [visual_files])
        .collect { visual_file -> "\"${visual_file.name}\"" }
        .join(', ')
    def prediction_file_names = (prediction_files instanceof List ? prediction_files : [prediction_files])
        .findAll { prediction_file -> prediction_file != null }
        .collect { prediction_file -> "\"${prediction_file.name}\"" }
        .join(', ')
"""
#!/usr/bin/env python3

import pandas as pd
import sys
import os
import shutil

def ordered_columns(frame, feature_list_path):
    if not feature_list_path or feature_list_path == 'null' or not os.path.exists(feature_list_path):
        return frame

    with open(feature_list_path) as handle:
        feature_order = [line.strip().split('\t')[0] for line in handle if line.strip()]

    ordered = [feature for feature in feature_order if feature in frame.columns]
    if ordered:
        return frame.loc[:, ordered]

    return frame

# Read DNAm and GEX transposed data (if provided)
def read_transposed(path):
    try:
        return pd.read_csv(path, sep=",", index_col=0)
    except FileNotFoundError:
        return pd.DataFrame()
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

dnam = read_transposed("${dnam_transposed}")
gex = read_transposed("${gex_transposed}")

# Read labels (classes)
try:
    labels_df = pd.read_csv("${classes_file}", sep="\t", index_col=0)
    labels = labels_df.iloc[:, 0]
except Exception:
    labels = pd.Series(dtype=object)

# Determine primary sample order: prefer DNAm, then GEX, then labels
if not dnam.empty:
    primary_samples = dnam.index.tolist()
elif not gex.empty:
    primary_samples = gex.index.tolist()
else:
    primary_samples = labels.index.tolist()

# Reindex labels to match primary sample order, but keep all samples even if labels are missing
labels = labels.reindex(primary_samples)

# Final sample order follows the data when available, otherwise the labels
final_samples = primary_samples if primary_samples else labels.index.tolist()

# Subset/reorder dataframes to final sample order where possible
if not dnam.empty:
    dnam = dnam.loc[[s for s in final_samples if s in dnam.index]]
    dnam = ordered_columns(dnam, "${params.dnam_probes_file ?: params.common_probes}")

if not gex.empty:
    gex = gex.loc[[s for s in final_samples if s in gex.index]]
    gex = ordered_columns(gex, "${params.gex_genes_file ?: ''}")

# Write labels as labels.csv with index renamed to "id"
labels.index.name = "id"
labels.to_csv("labels.csv", sep=",")

# Write DNAm features (preserve feature file order; samples now match labels)
if not dnam.empty:
    dnam.to_csv("features.dnam.csv", sep=",")
else:
    pd.DataFrame().to_csv("features.dnam.csv", sep=",")

# Write GEX features (preserve feature file order; samples now match labels)
if not gex.empty:
    gex.to_csv("features.gex.csv", sep=",")
else:
    pd.DataFrame().to_csv("features.gex.csv", sep=",")

# Write feature names for DNAM
if not dnam.empty:
    with open("feature_names.dnam.txt", "w") as f:
        f.write("\\n".join(dnam.columns.tolist()))
else:
    open("feature_names.dnam.txt", "w").close()

# Write feature names for GEX
if not gex.empty:
    with open("feature_names.gex.txt", "w") as f:
        f.write("\\n".join(gex.columns.tolist()))
else:
    open("feature_names.gex.txt", "w").close()

try:
    if '${gex_norm_factors_rds}' != 'null' and os.path.exists("${gex_norm_factors_rds}"):
        src = os.path.abspath("${gex_norm_factors_rds}")
        dst = os.path.abspath("gex_norm_factors.rds")
        if src != dst:
            shutil.copy2(src, dst)
    else:
        open('gex_norm_factors.rds', 'wb').close()
except Exception:
    open('gex_norm_factors.rds', 'wb').close()

# Write filtered class names if provided (one per line)
try:
    if '${gex_classes_filtered}' != 'null' and os.path.exists("${gex_classes_filtered}"):
        try:
            gex_classes_df = pd.read_csv("${gex_classes_filtered}", sep="\t", index_col=0)
            # If TSV has index=class and Count column, write index
            with open('class_names.gex.txt', 'w') as f:
                for cls in gex_classes_df.index.tolist():
                    f.write(f"{cls}\\n")
        except Exception:
            # Fallback: try reading as single-column file
            try:
                with open("${gex_classes_filtered}") as fin, open('class_names.gex.txt', 'w') as fout:
                    for line in fin:
                        parts = line.strip().split('\t')
                        if parts:
                            fout.write(parts[0] + '\\n')
            except Exception:
                open('class_names.gex.txt', 'w').close()
    else:
        open('class_names.gex.txt', 'w').close()
except Exception:
    open('class_names.gex.txt', 'w').close()

try:
    if '${dnam_classes_filtered}' != 'null' and os.path.exists("${dnam_classes_filtered}"):
        try:
            dnam_classes_df = pd.read_csv("${dnam_classes_filtered}", sep="\t", index_col=0)
            with open('class_names.dnam.txt', 'w') as f:
                for cls in dnam_classes_df.index.tolist():
                    f.write(f"{cls}\\n")
        except Exception:
            try:
                with open("${dnam_classes_filtered}") as fin, open('class_names.dnam.txt', 'w') as fout:
                    for line in fin:
                        parts = line.strip().split('\t')
                        if parts:
                            fout.write(parts[0] + '\\n')
            except Exception:
                open('class_names.dnam.txt', 'w').close()
    else:
        open('class_names.dnam.txt', 'w').close()
except Exception:
    open('class_names.dnam.txt', 'w').close()
os.makedirs("visuals", exist_ok=True)
for prediction_path in [${prediction_file_names}]:
    if os.path.exists(prediction_path):
        destination = os.path.basename(prediction_path)
        if os.path.abspath(prediction_path) != os.path.abspath(destination):
            shutil.copy2(prediction_path, destination)

visual_inputs = [${visual_file_names}]
for visual_path in visual_inputs:
    if not visual_path.lower().endswith((".png", ".svg")):
        continue

    destination = os.path.join("visuals", os.path.basename(visual_path))
    if not os.path.exists(destination):
        shutil.copy2(visual_path, destination)
        continue

    stem, ext = os.path.splitext(os.path.basename(visual_path))
    idx = 1
    while True:
        candidate = os.path.join("visuals", f"{stem}.{idx}{ext}")
        if not os.path.exists(candidate):
            shutil.copy2(visual_path, candidate)
            break
        idx += 1

# Create versions file
import pandas
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f'    python: "{sys.version.split()[0]}"\\n')
    f.write(f'    pandas: "{pandas.__version__}"\\n')
"""

    stub:
    """
    touch labels.csv
    touch features.dnam.csv
    touch features.gex.csv
    touch feature_names.dnam.txt
    touch feature_names.gex.txt
    touch gex_norm_factors.rds
    touch class_distribution.gex.filtered.svg
    touch class_distribution.dnam.filtered.svg
    touch classes.gex.filtered.tsv
    touch classes.dnam.filtered.tsv
    touch samplesheet.filtered.csv
    touch class_names.gex.txt
    touch class_names.dnam.txt
    mkdir -p visuals

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: "stub"
        pandas: "stub"
    END_VERSIONS
    """
}
