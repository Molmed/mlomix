#!/usr/bin/env python3

import argparse
from pathlib import Path

import joblib
import pandas as pd


def crc_threshold(calibration: pd.DataFrame, alpha: float) -> float:
    class_columns = [
        column
        for column in calibration.columns
        if column not in {"id", "true_label", "predicted_label"}
    ]
    probabilities = calibration[class_columns].apply(pd.to_numeric)
    true_probabilities = [
        probabilities.loc[index, true_label]
        for index, true_label in calibration["true_label"].items()
        if true_label in probabilities.columns
    ]
    if not true_probabilities:
        raise ValueError("Calibration data contains no usable true labels")
    n_calibration = len(true_probabilities)
    target_false_negative_rate = (
        (n_calibration + 1) / n_calibration * alpha
        - 1 / (n_calibration + 1)
    )
    target_false_negative_rate = max(0.0, min(1.0, target_false_negative_rate))
    return float(pd.Series(true_probabilities).quantile(target_false_negative_rate))


def predict(input_path: Path, model_path: Path, calibration_path: Path, alpha: float) -> pd.DataFrame:
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between 0 and 1")

    model = joblib.load(model_path)
    data = pd.read_csv(input_path, index_col=0)
    feature_names = list(getattr(model, "feature_names_in_", data.columns))
    missing = [feature for feature in feature_names if feature not in data.columns]
    if missing:
        raise ValueError(f"Input is missing {len(missing)} model features, including {missing[:3]}")

    probabilities = model.predict_proba(data.loc[:, feature_names])
    model_classes = [str(value) for value in model.classes_]
    predictions = pd.DataFrame(probabilities, index=data.index, columns=model_classes)
    calibration = pd.read_csv(calibration_path)
    calibration_classes = set(calibration["true_label"].dropna())
    missing_classes = sorted(set(model_classes) - calibration_classes)
    if missing_classes:
        raise ValueError(f"Calibration data is missing model classes: {missing_classes}")
    threshold = crc_threshold(calibration, alpha)

    result = pd.DataFrame(
        {
            "id": data.index,
            "predicted_label": [
                ";".join(predictions.columns[row >= threshold])
                for row in predictions.to_numpy()
            ],
        }
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--calibration", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()
    predict(args.input, args.model, args.calibration, args.alpha).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
