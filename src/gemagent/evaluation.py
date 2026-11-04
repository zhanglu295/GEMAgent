"""Final evaluation utilities.

Truth labels are used here only, after model construction. Runtime expansion and
prune modules never import this module.
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd


def _safe_div(num, denom):
    return float(num) / float(denom) if denom else 0.0


def confusion_metrics(tp, fp, tn, fn):
    sensitivity = _safe_div(tp, tp + fn)
    specificity = _safe_div(tn, tn + fp)
    precision = _safe_div(tp, tp + fp)
    recall = sensitivity
    f1 = _safe_div(2 * precision * recall, precision + recall)
    beta2 = 4.0
    f2 = _safe_div((1 + beta2) * precision * recall, beta2 * precision + recall)
    fpr = _safe_div(fp, fp + tn)
    balanced_accuracy = (sensitivity + specificity) / 2.0
    denom = math.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = _safe_div(tp * tn - fp * fn, denom)
    return {
        "n": tp + fp + tn + fn,
        "TP": tp, "FP": fp, "TN": tn, "FN": fn,
        "accuracy": _safe_div(tp + tn, tp + fp + tn + fn),
        "sensitivity": sensitivity,
        "specificity": specificity,
        "precision": precision,
        "recall": recall,
        "F1": f1,
        "F2": f2,
        "FPR": fpr,
        "balanced_accuracy": balanced_accuracy,
        "MCC": mcc,
    }


def decision_to_predicted(frame, decision_column="decision", applied_column="applied"):
    result = frame.copy()
    decisions = result[decision_column].astype(str).str.upper()
    applied = result[applied_column].astype(str).str.lower().isin(("true", "1", "yes"))
    result["predicted"] = 1
    result.loc[(decisions == "BLOCK") & applied, "predicted"] = 0
    return result


def evaluate_frames(truth, calls, predicted_column="predicted"):
    truth = truth.copy()
    calls = calls.copy()
    if "exchange" in truth.columns and "exchange_id" not in truth.columns:
        truth = truth.rename(columns={"exchange": "exchange_id"})
    if "exchange" in calls.columns and "exchange_id" not in calls.columns:
        calls = calls.rename(columns={"exchange": "exchange_id"})
    if "observed" not in truth.columns:
        raise ValueError("truth table requires an observed column")
    if predicted_column not in calls.columns:
        raise ValueError(f"calls table requires {predicted_column}")
    keys = [c for c in ("organism", "exchange_id", "direction") if c in truth.columns and c in calls.columns]
    if "exchange_id" not in keys:
        raise ValueError("truth and calls require a shared exchange_id column")
    keep = keys + [predicted_column]
    merged = truth.merge(calls[keep].drop_duplicates(), on=keys, how="left")
    merged[predicted_column] = pd.to_numeric(merged[predicted_column], errors="coerce").fillna(0).astype(int)
    observed = pd.to_numeric(merged["observed"], errors="coerce").fillna(0).astype(int)
    pred = merged[predicted_column].astype(int)
    tp = int(((pred == 1) & (observed == 1)).sum())
    fp = int(((pred == 1) & (observed == 0)).sum())
    tn = int(((pred == 0) & (observed == 0)).sum())
    fn = int(((pred == 0) & (observed == 1)).sum())
    return confusion_metrics(tp, fp, tn, fn), merged


def evaluate_files(truth_path, calls_path, predicted_column="predicted"):
    truth = pd.read_csv(Path(truth_path), sep="\t", dtype={"exchange_id": str})
    calls = pd.read_csv(Path(calls_path), sep="\t", dtype={"exchange_id": str})
    return evaluate_frames(truth, calls, predicted_column=predicted_column)
