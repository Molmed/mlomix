#!/usr/bin/env python3

import argparse
import ast
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from upsetplot import plot


def scalarize_text_positions(fig):
    for axes in fig.axes:
        for text in axes.texts:
            x_position, y_position = text.get_position()
            if isinstance(x_position, np.ndarray) and x_position.size == 1:
                x_position = x_position.item()
            if isinstance(y_position, np.ndarray) and y_position.size == 1:
                y_position = y_position.item()
            text.set_position((x_position, y_position))


def normalize_label(value):
    if pd.isna(value):
        return "EMPTY SET"
    value = str(value).strip()
    return "EMPTY SET" if not value or value.lower() == "nan" else value


def parse_prediction_classes(value):
    if isinstance(value, (list, tuple, set, frozenset)):
        return [normalize_label(item) for item in value]
    value = str(value).strip()
    try:
        parsed = ast.literal_eval(value)
    except (ValueError, SyntaxError):
        parsed = value.split(";")
    if isinstance(parsed, (list, tuple, set, frozenset)):
        return [normalize_label(item) for item in parsed]
    return [normalize_label(parsed)]


def make_upset(prediction_sets, output):
    class_names = sorted(
        {label for values in prediction_sets for label in values},
        key=str,
    )
    upset_data = pd.DataFrame(
        [
            {class_name: class_name in values for class_name in class_names}
            for values in prediction_sets
        ],
        dtype=bool,
    )
    if class_names:
        upset_data = upset_data.loc[:, upset_data.any(axis=0)]
        upset_data.set_index(upset_data.columns.tolist(), inplace=True)
    fig = plt.figure()
    plot(
        upset_data,
        fig=fig,
        element_size=30,
        sort_by="cardinality",
        facecolor="#9768ac",
        show_counts="%d",
        show_percentages=False,
    )
    scalarize_text_positions(fig)
    plt.tight_layout()
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)


def make_sankey(data, colors_path, output):
    palette_df = pd.read_csv(
        colors_path,
        sep="\t",
        header=None,
        names=["class", "color"],
        index_col=0,
    )
    flow_counts = defaultdict(int)
    for _, row in data.iterrows():
        true_label = normalize_label(row["true_label"])
        for predicted_label in parse_prediction_classes(row["predicted_label"]):
            flow_counts[(true_label, predicted_label)] += 1

    true_labels = sorted({key[0] for key in flow_counts})
    predicted_labels = sorted({key[1] for key in flow_counts})
    suffix = " "
    display_predictions = [f"{label}{suffix}" for label in predicted_labels]
    all_labels = true_labels + display_predictions
    label_index = {label: index for index, label in enumerate(all_labels)}

    sources, targets, values = [], [], []
    for (true_label, predicted_label), count in flow_counts.items():
        sources.append(label_index[true_label])
        targets.append(label_index[f"{predicted_label}{suffix}"])
        values.append(count)

    palette = {
        str(class_name).strip(): str(color).strip()
        for class_name, color in palette_df["color"].items()
    }
    node_colors = [palette.get(label.rstrip(), "#BDBDBD") for label in all_labels]

    def rgba(color, alpha):
        color = color.lstrip("#")
        return f"rgba({int(color[0:2], 16)},{int(color[2:4], 16)},{int(color[4:6], 16)},{alpha})"

    figure = go.Figure(
        go.Sankey(
            arrangement="snap",
            node=dict(
                pad=15,
                thickness=15,
                line=dict(color="black", width=0.5),
                label=all_labels,
                color=node_colors,
            ),
            link=dict(
                source=sources,
                target=targets,
                value=values,
                color=[rgba(node_colors[source], 0.5) for source in sources],
            ),
        )
    )
    figure.update_layout(
        title_text="Conformal prediction sets",
        height=500 + 40 * max(len(true_labels), len(predicted_labels)),
        margin=dict(l=20, r=20, t=60, b=20),
        font=dict(family="Arial, sans-serif", size=16, color="black"),
    )
    figure.write_html(output, include_plotlyjs=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--labels", type=Path, required=True)
    parser.add_argument("--colors", type=Path, required=True)
    parser.add_argument("--modality", required=True)
    args = parser.parse_args()

    predictions = pd.read_csv(args.predictions)
    labels = pd.read_csv(args.labels, sep="\t", index_col=0)
    labels = labels.iloc[:, 0].rename("true_label")
    combined = predictions.set_index("id").join(labels, how="left").reset_index()
    prediction_sets = [
        set(parse_prediction_classes(value))
        for value in combined["predicted_label"]
    ]
    make_upset(prediction_sets, f"prediction_sets.{args.modality}.upset.png")
    make_sankey(
        combined,
        args.colors,
        f"prediction_sets.{args.modality}.sankey.html",
    )


if __name__ == "__main__":
    main()
