# -*- coding: utf-8 -*-
"""Script 7: OSBs recall analysis

Purpose
-------
Analyze independent repeat-identification and observer agreement.

Method overview
---------------
1. Restrict the independent repeat-identification table to the final superbubble sample.
2. Normalize the binary observer decisions and calculate object recall, pairwise
   agreement and multi-observer agreement statistics.
3. Write the object-level and aggregate results and their summary figure; the
   observations here are reader decisions, not simulated trials.

Main inputs
-----------
- ../data/third_party_recall_test.csv
- ../results/superbubble_final_fit_parameters.csv

Main outputs
------------
- ../results/figures/7_independent_repeat_recall.png; independent-identification recall figure.

Figure/table role
-----------------
../results/figures/7_independent_repeat_recall.png; independent-identification recall figure.

Runtime and data notes
----------------------
Reference runtime: 0.95 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

External-data caveat
--------------------
This script uses packaged intermediate or final inputs unless the inputs section explicitly names an external cache or survey product.

Reading the code
----------------
Start with the path and scientific settings below, then follow main() at
the end of the file. The method overview above describes the order of the
analysis steps; the input/output lists identify the upstream data and
products written by this stage. Relative data paths are resolved from code/.
"""
from __future__ import annotations
import time

import argparse
import itertools
import json
import os
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import pandas as pd


# Path configuration; paths are relative to code/ and point directly into ../data or ../results/.
HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "7_independent_repeat_recall"
FINAL_FIG_DIR = Path("..") / "results" / "figures"

DEFAULT_REPEAT_CSV = DATA_DIR / "third_party_recall_test.csv"
DEFAULT_FINAL_SAMPLE_CSV = Path("..") / "results" / "superbubble_final_fit_parameters.csv"

PRIMARY_READER_COLS = ["wt", "hbw", "ggr"]
INDEPENDENT_READER_COLS = ["lfp", "wj", "sly", "zzh", "zjm", "lgh"]
ALL_READER_COLS = PRIMARY_READER_COLS + INDEPENDENT_READER_COLS
MANUAL_NEW_CANDIDATE_IDS = [34, 35]


def ensure_existing_file(path, label):
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"{label}  does not exist: {p}")
    return p


def load_final_sample_ids(path):
    df = pd.read_csv(path)
    if "id" not in df.columns:
        raise ValueError('Invalid input or missing required data.')
    if "selected_ge_4" in df.columns:
        df = df[df["selected_ge_4"].astype(bool)]
    return [int(v) for v in df["id"].tolist()]


def normalize_tf(value):
    return str(value).strip().upper() == "T"


def load_repeat_table(path, reader_cols):
    df = pd.read_csv(path)
    if "ID" not in df.columns:
        raise ValueError('Invalid input or missing required data.')
    missing = [col for col in reader_cols if col not in df.columns]
    if missing:
        raise ValueError('Invalid input or missing required data.')
    for col in reader_cols:
        df[col] = df[col].map(normalize_tf)
    df["ID"] = df["ID"].astype(int)
    return df


def add_object_metrics(df, reader_cols):
    out = df.copy()
    out["n_readers"] = len(reader_cols)
    out["n_detected"] = out[reader_cols].sum(axis=1).astype(int)
    out["recall_fraction"] = out["n_detected"] / len(reader_cols)
    out["recall_percent"] = 100.0 * out["recall_fraction"]
    out["unanimous_detected"] = out["n_detected"] == len(reader_cols)
    out["majority_detected"] = out["n_detected"] >= (len(reader_cols) // 2 + 1)
    out["at_least_one_detected"] = out["n_detected"] >= 1
    return out


def pairwise_metrics(main_df, reader_cols):
    rows = []
    for a, b in itertools.combinations(reader_cols, 2):
        x = main_df[a].to_numpy(dtype=bool)
        y = main_df[b].to_numpy(dtype=bool)
        both_t = int(np.sum(x & y))
        either_t = int(np.sum(x | y))
        rows.append(
            {
                "reader_a": a,
                "reader_b": b,
                "agreement_fraction": float(np.mean(x == y)),
                "jaccard_detected_fraction": float(both_t / either_t) if either_t else float("nan"),
                "both_detected_count": both_t,
                "either_detected_count": either_t,
            }
        )
    return pd.DataFrame(rows)


def fleiss_kappa_binary(object_df, reader_cols):
    n = len(reader_cols)
    counts_t = object_df["n_detected"].to_numpy(dtype=float)
    counts = np.column_stack([n - counts_t, counts_t])
    if len(counts) == 0:
        return float("nan")
    p = counts.sum(axis=0) / (len(counts) * n)
    p_i = ((counts**2).sum(axis=1) - n) / (n * (n - 1))
    p_bar = float(np.mean(p_i))
    p_e = float(np.sum(p**2))
    if abs(1.0 - p_e) < 1e-12:
        return float("nan")
    return float((p_bar - p_e) / (1.0 - p_e))


def build_summary(main_objects, nonretained_candidates, observer_summary, pairwise_df, reader_cols):
    total_detected = int(main_objects["n_detected"].sum())
    total_trials = int(len(main_objects) * len(reader_cols))
    recall_values = main_objects["recall_fraction"].to_numpy(dtype=float)
    summary = {
        "reader_columns_used": reader_cols,
        "primary_reader_columns_shown": PRIMARY_READER_COLS,
        "n_retained_targets": int(len(main_objects)),
        "n_independent_readers": int(len(reader_cols)),
        "total_reader_target_decisions": total_trials,
        "total_detected_decisions": total_detected,
        "overall_individual_recall_fraction": float(total_detected / total_trials),
        "mean_object_recall_fraction": float(np.mean(recall_values)),
        "median_object_recall_fraction": float(np.median(recall_values)),
        "min_object_recall_fraction": float(np.min(recall_values)),
        "unanimous_object_fraction": float(np.mean(main_objects["unanimous_detected"])),
        "object_fraction_detected_by_at_least_5_of_6": float(np.mean(main_objects["n_detected"] >= 5)),
        "object_fraction_detected_by_majority": float(np.mean(main_objects["majority_detected"])),
        "object_fraction_detected_by_at_least_one": float(np.mean(main_objects["at_least_one_detected"])),
        "mean_reader_recall_fraction": float(observer_summary["recall_fraction"].mean()),
        "min_reader_recall_fraction": float(observer_summary["recall_fraction"].min()),
        "max_reader_recall_fraction": float(observer_summary["recall_fraction"].max()),
        "mean_pairwise_agreement_fraction": float(pairwise_df["agreement_fraction"].mean()),
        "min_pairwise_agreement_fraction": float(pairwise_df["agreement_fraction"].min()),
        "mean_pairwise_detected_jaccard_fraction": float(pairwise_df["jaccard_detected_fraction"].mean()),
        "fleiss_kappa_binary": fleiss_kappa_binary(main_objects, reader_cols),
        "manual_new_candidate_ids_not_retained": MANUAL_NEW_CANDIDATE_IDS,
        "nonretained_candidate_ids": [int(v) for v in nonretained_candidates["ID"].tolist()],
        "nonretained_candidate_support": {
            f"SB{int(row.ID)}": {
                "n_detected": int(row.n_detected),
                "support_fraction": float(row.recall_fraction),
            }
            for row in nonretained_candidates.itertuples(index=False)
        },
    }
    return summary


def plot_summary(main_objects, nonretained_candidates, observer_summary, summary, reader_cols, output_path):
    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.titlesize": 10,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 7.5,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    main_objects = main_objects.sort_values("ID").reset_index(drop=True)
    all_reader_cols = PRIMARY_READER_COLS + reader_cols
    raw_matrix = main_objects[all_reader_cols].to_numpy(dtype=int)
    matrix = raw_matrix.copy()
    matrix[:, :len(PRIMARY_READER_COLS)] *= 2
    y = np.arange(len(main_objects))
    fig = plt.figure(figsize=(7.8, 7.2), constrained_layout=True)
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[5.35, 1.35],
        width_ratios=[1.35, 1.0],
        hspace=0.08,
        wspace=0.12,
    )
    ax_heat = fig.add_subplot(gs[0, 0])
    ax_obj = fig.add_subplot(gs[0, 1], sharey=ax_heat)
    ax_reader = fig.add_subplot(gs[1, 0])
    ax_new = fig.add_subplot(gs[1, 1])

    cmap = ListedColormap(["#FFFFFF", "#2A9D8F", "#4E79A7"])
    ax_heat.imshow(matrix, aspect="auto", cmap=cmap, vmin=0, vmax=2, interpolation="nearest")
    ax_heat.set_xticks(np.arange(len(all_reader_cols)), all_reader_cols, rotation=45, ha="right")
    for col_i, tick_label in enumerate(ax_heat.get_xticklabels()):
        tick_label.set_color("#4E79A7" if col_i < len(PRIMARY_READER_COLS) else "#2A9D8F")
    ax_heat.set_yticks(y, [f"SB{int(v)}" for v in main_objects["ID"]])
    ax_heat.set_title("Independent recognition matrix", pad=6)
    ax_heat.set_xlabel("Reader")
    ax_heat.set_ylabel("Final sample")
    ax_heat.set_xticks(np.arange(-0.5, len(all_reader_cols), 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, len(main_objects), 1), minor=True)
    ax_heat.grid(which="minor", color="0.82", linewidth=0.45)
    ax_heat.tick_params(which="minor", bottom=False, left=False)
    for row_i in range(matrix.shape[0]):
        for col_i in range(matrix.shape[1]):
            if matrix[row_i, col_i] == 0:
                ax_heat.text(col_i, row_i, "x", ha="center", va="center", color="#C44E52", fontsize=7.5, fontweight="bold")
    ax_heat.axvline(len(PRIMARY_READER_COLS) - 0.5, color="0.15", linewidth=1.0)

    recall_pct = 100.0 * main_objects["recall_fraction"].to_numpy(dtype=float)
    colors = np.where(recall_pct >= 100.0, "#2A9D8F", np.where(recall_pct >= 80.0, "#E9C46A", "#E76F51"))
    ax_obj.barh(y, recall_pct, color=colors, edgecolor="0.15", linewidth=0.45)
    ax_obj.axvline(100.0 * summary["overall_individual_recall_fraction"], color="black", linewidth=0.95)
    ax_obj.set_xlim(0, 105)
    ax_obj.set_xticks([0, 50, 100])
    ax_obj.set_xlabel("Object recall (%)")
    ax_obj.set_title("Recall per target", pad=6)
    ax_obj.tick_params(axis="y", labelleft=False, left=False)
    ax_obj.grid(True, axis="x", color="0.88", linewidth=0.55)
    for yi, value, count in zip(y, recall_pct, main_objects["n_detected"]):
        ax_obj.text(
            -0.035,
            yi,
            f"{int(count)}/6",
            transform=ax_obj.get_yaxis_transform(),
            va="center",
            ha="right",
            fontsize=6.8,
            clip_on=False,
        )
    ax_obj.text(
        100.0 * summary["overall_individual_recall_fraction"] - 1.0,
        -0.75,
        f"overall {100.0 * summary['overall_individual_recall_fraction']:.1f}%",
        ha="right",
        va="bottom",
        rotation=90,
        fontsize=6.7,
        color="0.15",
    )

    reader_pct = 100.0 * observer_summary["recall_fraction"].to_numpy(dtype=float)
    x_reader = np.arange(len(observer_summary))
    reader_colors = [
        "#4E79A7" if reader in PRIMARY_READER_COLS else "#2A9D8F"
        for reader in observer_summary["reader"].tolist()
    ]
    ax_reader.bar(x_reader, reader_pct, color=reader_colors, edgecolor="0.15", linewidth=0.45)
    ax_reader.axhline(100.0 * summary["overall_individual_recall_fraction"], color="black", linewidth=0.85)
    ax_reader.set_ylim(0, 105)
    ax_reader.set_yticks([0, 50, 100])
    ax_reader.set_xticks(x_reader, observer_summary["reader"].tolist(), rotation=45, ha="right")
    ax_reader.set_ylabel("Reader support (%)")
    ax_reader.set_title("Recall per reader", pad=5)
    ax_reader.grid(True, axis="y", color="0.88", linewidth=0.55)
    for xi, value in zip(x_reader, reader_pct):
        ax_reader.text(xi, min(value + 3, 102), f"{value:.0f}", ha="center", va="bottom", fontsize=6.5)

    nonretained_sorted = nonretained_candidates.sort_values("ID").reset_index(drop=True)
    nonretained_pct = 100.0 * nonretained_sorted["recall_fraction"].to_numpy(dtype=float)
    x_new = np.arange(len(nonretained_sorted))
    bar_colors = np.where(
        nonretained_sorted["ID"].isin(MANUAL_NEW_CANDIDATE_IDS),
        "#B8B8B8",
        "#D8C7A2",
    )
    ax_new.bar(x_new, nonretained_pct, color=bar_colors, edgecolor="0.15", linewidth=0.45)
    ax_new.set_ylim(0, 105)
    ax_new.set_yticks([0, 50, 100])
    ax_new.set_xticks(x_new, [f"SB{int(v)}" for v in nonretained_sorted["ID"]], rotation=45, ha="right")
    ax_new.set_ylabel("Support (%)")
    ax_new.set_title("Objects not in final sample", pad=5)
    ax_new.grid(True, axis="y", color="0.88", linewidth=0.55)
    for xi, value, count in zip(x_new, nonretained_pct, nonretained_sorted["n_detected"]):
        ax_new.text(xi, min(value + 4, 101), f"{int(count)}/6", ha="center", va="bottom", fontsize=6.3)

    for ax in [ax_heat, ax_obj, ax_reader, ax_new]:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    ax_heat.spines["right"].set_visible(True)
    ax_heat.spines["top"].set_visible(True)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=260)
    plt.close(fig)


def build_parser():
    """Build the command-line interface for Script 7."""
    parser = argparse.ArgumentParser(description='Run Script 7: OSBs recall analysis.')
    parser.add_argument("--repeat-csv", default=str(DEFAULT_REPEAT_CSV), help='Input file used by this stage.')
    parser.add_argument("--final-sample-csv", default=str(DEFAULT_FINAL_SAMPLE_CSV), help='Input file used by this stage.')
    return parser


def main():
    """Run Script 7 from validated inputs to the documented outputs."""
    args = build_parser().parse_args()
    ensure_existing_file(args.repeat_csv, "repeat-reader candidate table")
    ensure_existing_file(args.final_sample_csv, "final publication sample table")
    final_sample_ids = load_final_sample_ids(args.final_sample_csv)
    repeat_df = load_repeat_table(args.repeat_csv, ALL_READER_COLS)

    main_objects = add_object_metrics(repeat_df[repeat_df["ID"].isin(final_sample_ids)], INDEPENDENT_READER_COLS)
    if len(main_objects) != len(final_sample_ids):
        missing = sorted(set(final_sample_ids) - set(main_objects["ID"].tolist()))
        raise RuntimeError('Runtime failure in this pipeline stage.')
    nonretained_candidates = add_object_metrics(repeat_df[~repeat_df["ID"].isin(final_sample_ids)], INDEPENDENT_READER_COLS)

    observer_rows = []
    for reader in ALL_READER_COLS:
        count = int(main_objects[reader].sum())
        observer_rows.append(
            {
                "reader": reader,
                "detected_count": count,
                "total_retained_targets": int(len(main_objects)),
                "recall_fraction": float(count / len(main_objects)),
                "recall_percent": float(100.0 * count / len(main_objects)),
            }
        )
    observer_summary = pd.DataFrame(observer_rows)
    pairwise_df = pairwise_metrics(main_objects, INDEPENDENT_READER_COLS)
    summary = build_summary(main_objects, nonretained_candidates, observer_summary, pairwise_df, INDEPENDENT_READER_COLS)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)
    object_csv = OUT_DIR / "7_recall_by_target_details.csv"
    observer_csv = OUT_DIR / "7_recall_by_identifier.csv"
    pairwise_csv = OUT_DIR / "7_pairwise_agreement.csv"
    summary_json = OUT_DIR / "7_overall_statistics.json"
    figure_path = FINAL_FIG_DIR / "7_independent_repeat_recall.png"

    main_objects.to_csv(object_csv, index=False)
    observer_summary.to_csv(observer_csv, index=False)
    pairwise_df.to_csv(pairwise_csv, index=False)
    summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    plot_summary(main_objects, nonretained_candidates, observer_summary, summary, INDEPENDENT_READER_COLS, figure_path)

    print("Repeat-reader recall analysis complete.")
    print(f"Final-sample objects: {len(main_objects)}")
    print(f"Non-retained candidates: {len(nonretained_candidates)}")
    print(f"Independent readers: {len(INDEPENDENT_READER_COLS)}")
    print(f"Figure: {figure_path}")
    print(f"Fleiss kappa = {summary['fleiss_kappa_binary']:.3f}")
    print(f"recall_by_target_details = {object_csv}")
    print(f"recall_by_identifier   = {observer_csv}")
    print(f"pairwise_agreement     = {pairwise_csv}")
    print(f"overall_statistics json  = {summary_json}")
    print("Script 7 complete.")

def _format_runtime(seconds: float) -> str:
    """Return a compact human-readable wall-clock runtime string."""
    seconds = max(0.0, float(seconds))
    hours, remainder = divmod(seconds, 3600.0)
    minutes, seconds = divmod(remainder, 60.0)
    if hours >= 1.0:
        return f"{int(hours)} h {int(minutes):02d} min {seconds:05.2f} s"
    if minutes >= 1.0:
        return f"{int(minutes)} min {seconds:05.2f} s"
    return f"{seconds:.2f} s"


def _run_with_timing(entrypoint, script_file: str) -> None:
    """Run a script entry point and print its total wall-clock runtime."""
    script_name = Path(script_file).name
    start = time.perf_counter()
    try:
        entrypoint()
    except SystemExit as exc:
        elapsed = time.perf_counter() - start
        if exc.code in (None, 0):
            print(f"[runtime] {script_name} completed in: {_format_runtime(elapsed)} ({elapsed:.2f} s)", flush=True)
        else:
            print(f"[runtime] {script_name} elapsed time before failure: {_format_runtime(elapsed)} ({elapsed:.2f} s)", flush=True)
        raise
    except BaseException:
        elapsed = time.perf_counter() - start
        print(f"[runtime] {script_name} elapsed time before failure: {_format_runtime(elapsed)} ({elapsed:.2f} s)", flush=True)
        raise
    elapsed = time.perf_counter() - start
    print(f"[runtime] {script_name} completed in: {_format_runtime(elapsed)} ({elapsed:.2f} s)", flush=True)


if __name__ == "__main__":
    _run_with_timing(main, __file__)
