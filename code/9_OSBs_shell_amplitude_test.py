"""Script 9: OSBs shell amplitude test

Purpose
-------
Compare shell jump amplitudes with extinction-map and stellar-sample uncertainties; automatically detects an installed dustmaps3d cache when needed.

Method overview
---------------
1. Load the final shell geometry, the packaged extinction-error map and the installed
   dustmaps3d component data.
2. Scan the lines of sight intersecting each shell, associate fitted extinction
   components and attach map and stellar-sample error estimates.
3. Summarize the shell amplitudes and their three-sigma comparisons, then export the
   target-level tables and figure.

Main inputs
-----------
- ../results/superbubble_final_fit_parameters.csv
- ../data/extinction_error.fits
- dustmaps3d local cache data_v3.fits

Main outputs
------------
- ../results/figures/9_shell_amplitude_vs_extinction_error_3sigma.png; shell-significance figure.

Figure/table role
-----------------
../results/figures/9_shell_amplitude_vs_extinction_error_3sigma.png; shell-significance figure.

Runtime and data notes
----------------------
Reference runtime: 35.22 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

External-data caveat
--------------------
A full rerun requires a local dustmaps3d installation and cache; the script auto-detects that cache when installed.

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
import importlib.util
import math
import os
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.io import fits
from astropy_healpix import HEALPix

# ----------------------------------------------------------------------------
# Path configuration; paths are relative to code/ and point directly into ../data or ../results/.
# ----------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "9_shell_amplitude_error_analysis"

GEOM_TABLE = FINAL_DIR / "superbubble_final_fit_parameters.csv"
ERROR_FITS = DATA_DIR / "extinction_error.fits"

OUT_SUMMARY = OUT_DIR / "9_shell_amplitude_error_statistics.csv"
OUT_DETAIL = OUT_DIR / "9_shell component details.csv"
OUT_FIG = FINAL_FIG_DIR / "9_shell_amplitude_vs_extinction_error_3sigma.png"

# ----------------------------------------------------------------------------
# Geometry and shell-mask conventions used by this analysis stage.
# ----------------------------------------------------------------------------
NSIDE = 1024
SHELL_U_MIN = 0.9
SHELL_U_MAX = 1.2
CYLINDER_HALF_HEIGHT_KPC = 0.025
AMP_ERROR_RATIO_THRESHOLD = 3.0
CHUNK_SIZE = 500_000

HEALPIX = HEALPix(nside=NSIDE, order="ring")


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
def resolve_data_fits(explicit: Path | None) -> Path:
    """Helper for the documented workflow; scientific assumptions are described in the module header."""
    if explicit is not None:
        if not explicit.exists():
            raise FileNotFoundError('Missing input file.')
        return explicit

    if importlib.util.find_spec("dustmaps3d") is None:
        raise ImportError(
            "dustmaps3d is required when --data-fits is not supplied; install it "
            "or provide a local package with --data-fits."
        )

    from platformdirs import user_data_dir

    cached = Path(user_data_dir("dustmaps3d")) / "data_v3.fits"
    if cached.exists():
        return cached

    try:
        from dustmaps3d import dustmaps3d

        dustmaps3d(0.0, 0.0, 1.0)
    except Exception as exc:
        raise FileNotFoundError(
            "dustmaps3d is installed, but its cached data_v3.fits file could not "
            "be found or initialized automatically."
        ) from exc

    if cached.exists():
        return cached

    raise FileNotFoundError(
        "dustmaps3d did not expose a cached data_v3.fits file after automatic "
        "initialization."
    )


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
def load_bubbles(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path)
    raw.columns = [c.lstrip("﻿") for c in raw.columns]
    df = pd.DataFrame()
    df["id"] = raw["id"].astype(int)
    df["mark"] = raw["mark"].astype(int)
    df["shape"] = raw["shape"].astype(str)
    df["center_x_kpc"] = raw["center_x_kpc"].astype(float)
    df["center_y_kpc"] = raw["center_y_kpc"].astype(float)
    df["center_z_kpc"] = raw["center_z_kpc"].astype(float)
    df["a_kpc"] = raw["a_radius_kpc"].astype(float)
    df["b_kpc"] = raw["b_radius_kpc"].astype(float)
    df["c_kpc"] = raw["c_radius_kpc"].astype(float)
    df["pa_deg"] = raw["angle_deg"].astype(float)
    df["xy_plane_z_kpc"] = raw["xy_plane_z_kpc"].astype(float)
    df["xy_plane_a_kpc"] = raw["xy_plane_a_kpc"].astype(float)
    df["xy_plane_b_kpc"] = raw["xy_plane_b_kpc"].astype(float)
    df["fit_std_pc"] = raw["fit_std_pc"].astype(float)
    df["source_type"] = np.where(
        raw["final_parameter_estimator"].astype(str) == "manual_override", "manual", "auto"
    )
    return df.sort_values("id").reset_index(drop=True)


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
def sigmoid_parameters(bubble, distance, span, cum_ebv):
    center = distance + 2.0 * span + bubble
    slope = 5.0 / span
    exponent = np.clip(5.0 * center / span, -700.0, 700.0)
    amplitude = cum_ebv * (1.0 / np.exp(exponent) + 1.0)
    return center, slope, amplitude


def rotate_to_local(dx, dy, pa_deg):
    """Rotate coordinates into the local frame used by the shell or projection calculation."""
    theta = math.radians(pa_deg)
    cos_t, sin_t = math.cos(theta), math.sin(theta)
    return dx * cos_t + dy * sin_t, -dx * sin_t + dy * cos_t


# ----------------------------------------------------------------------------
# Geometry and shell-mask conventions used by this analysis stage.
# ----------------------------------------------------------------------------
def shell_mask_cylinder(x, y, z, b):
    dx = x - b["center_x_kpc"]
    dy = y - b["center_y_kpc"]
    dz = z - b["center_z_kpc"]
    lx, ly = rotate_to_local(dx, dy, b["pa_deg"])
    a = max(b["a_kpc"], 1e-9)
    bb = max(b["b_kpc"], 1e-9)
    q = (lx / a) ** 2 + (ly / bb) ** 2
    mask = (np.abs(dz) <= CYLINDER_HALF_HEIGHT_KPC) & (q >= SHELL_U_MIN ** 2) & (q <= SHELL_U_MAX ** 2)
    return mask


def shell_mask_ellipsoid_cap(x, y, z, b):
    cz, c_full, xy_z = b["center_z_kpc"], b["c_kpc"], b["xy_plane_z_kpc"]
    mark = int(b["mark"])
    z_apex = cz - c_full if mark == 1 else cz + c_full
    c_cap = max(abs(xy_z - z_apex), 1e-9)
    cap_sign = -1.0 if mark == 1 else 1.0
    A = max(b["xy_plane_a_kpc"], 1e-9)
    B = max(b["xy_plane_b_kpc"], 1e-9)
    dx = x - b["center_x_kpc"]
    dy = y - b["center_y_kpc"]
    Z = z - xy_z
    lx, ly = rotate_to_local(dx, dy, b["pa_deg"])
    q = (lx / A) ** 2 + (ly / B) ** 2 + (Z / c_cap) ** 2
    cap_side = (cap_sign * Z) >= 0.0
    return cap_side & (q >= SHELL_U_MIN ** 2) & (q <= SHELL_U_MAX ** 2)


def shell_mask(x, y, z, b):
    if b["shape"] == "cylinder":
        return shell_mask_cylinder(x, y, z, b)
    return shell_mask_ellipsoid_cap(x, y, z, b)


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
def assign_distance_binned_min(distance, full, v02, v14, v2max):
    out = np.full(distance.shape, np.nan, dtype=float)
    m = distance < 1.0
    if np.any(m):
        out[m] = np.nanmin(np.vstack([full[m], v02[m]]), axis=0)
    m = (distance >= 1.0) & (distance < 2.0)
    if np.any(m):
        out[m] = np.nanmin(np.vstack([full[m], v02[m], v14[m]]), axis=0)
    m = (distance >= 2.0) & (distance < 4.0)
    if np.any(m):
        primary = np.nanmin(np.vstack([v14[m], v2max[m]]), axis=0)
        out[m] = np.where(np.isnan(primary), full[m], primary)
    m = distance >= 4.0
    if np.any(m):
        out[m] = np.where(np.isnan(v2max[m]), full[m], v2max[m])
    return out


# ----------------------------------------------------------------------------
# Geometry and shell-mask conventions used by this analysis stage.
# ----------------------------------------------------------------------------
def pixel_unit_vectors(pixel_ids):
    lon, lat = HEALPIX.healpix_to_lonlat(pixel_ids)
    l = lon.to_value(u.rad)
    bb = lat.to_value(u.rad)
    cos_b = np.cos(bb)
    return cos_b * np.cos(l), cos_b * np.sin(l), np.sin(bb), lon.to_value(u.deg), lat.to_value(u.deg)


def scan_shell_components(bubbles, data_fits, chunk_size, allow_beyond_max, quiet):
    records: list[dict] = []
    bubble_dicts = [row.to_dict() for _, row in bubbles.iterrows()]
    with fits.open(data_fits, memmap=True) as hdul:
        table = hdul[1].data
        n_rows = len(table)
        for start in range(0, n_rows, chunk_size):
            end = min(start + chunk_size, n_rows)
            pixel_ids = np.arange(start, end, dtype=np.int64)
            ux, uy, uz, l_deg, b_deg = pixel_unit_vectors(pixel_ids)
            bubble_col = np.asarray(table["bubble"][start:end], dtype=float)
            max_distance = np.asarray(table["max_distance"][start:end], dtype=float)
            for comp in range(1, 5):
                distance = np.asarray(table[f"distance_{comp}"][start:end], dtype=float)
                span = np.asarray(table[f"span_{comp}"][start:end], dtype=float)
                cum_ebv = np.asarray(table[f"Cum_EBV_{comp}"][start:end], dtype=float)
                valid = (
                    np.isfinite(distance) & np.isfinite(span) & np.isfinite(cum_ebv)
                    & (span > 0.0) & (cum_ebv != 0.0)
                )
                if not np.any(valid):
                    continue
                center, slope, amplitude = sigmoid_parameters(bubble_col, distance, span, cum_ebv)
                valid &= np.isfinite(center) & np.isfinite(slope) & np.isfinite(amplitude) & (center > 0.0)
                if not allow_beyond_max:
                    valid &= center <= max_distance
                if not np.any(valid):
                    continue
                idx = np.flatnonzero(valid)
                x = ux[idx] * center[idx]
                y = uy[idx] * center[idx]
                z = uz[idx] * center[idx]
                for b in bubble_dicts:
                    mask = shell_mask(x, y, z, b)
                    if not np.any(mask):
                        continue
                    glob = idx[np.flatnonzero(mask)]
                    for g in glob:
                        records.append(
                            {
                                "id": b["id"],
                                "component": comp,
                                "healpix_pixel": int(pixel_ids[g]),
                                "l_deg": float(l_deg[g]),
                                "b_deg": float(b_deg[g]),
                                "d_kpc": float(center[g]),
                                "sigmoid_amplitude_mag": float(amplitude[g]),
                                "shape": b["shape"],
                                "mark": int(b["mark"]),
                                "source_type": b["source_type"],
                            }
                        )
            if not quiet:
                print(
                    f"Scanned dust rows {start}-{end} of {n_rows}; shell-component hits so far={len(records)}",
                    flush=True,
                )
    return pd.DataFrame(records)


def attach_error_and_sigma(detail, data_fits, error_fits):
    if detail.empty:
        detail["extinction_error_mag"] = []
        detail["stellar_sample_sigma_mag"] = []
        return detail
    pixel = detail["healpix_pixel"].to_numpy(np.int64)
    distance = detail["d_kpc"].to_numpy(float)
    with fits.open(error_fits, memmap=True) as h:
        t = h[1].data
        err = assign_distance_binned_min(
            distance,
            np.asarray(t["extinction_error"][pixel], float),
            np.asarray(t["extinction_error02"][pixel], float),
            np.asarray(t["extinction_error14"][pixel], float),
            np.asarray(t["extinction_error2max"][pixel], float),
        )
    with fits.open(data_fits, memmap=True) as h:
        t = h[1].data
        sig = assign_distance_binned_min(
            distance,
            np.asarray(t["sigma"][pixel], float),
            np.asarray(t["sigma_0_2"][pixel], float),
            np.asarray(t["sigma_1_4"][pixel], float),
            np.asarray(t["sigma_2_max"][pixel], float),
        )
    detail = detail.copy()
    detail["extinction_error_mag"] = err
    detail["stellar_sample_sigma_mag"] = sig
    return detail


# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
def _mean_positive(series):
    v = pd.to_numeric(series, errors="coerce").to_numpy(float)
    v = v[np.isfinite(v) & (v > 0.0)]
    return float(np.mean(v)) if v.size else np.nan


def summarize(bubbles, detail):
    rows = []
    for _, b in bubbles.iterrows():
        bid = int(b["id"])
        sub = detail[detail["id"] == bid] if not detail.empty else detail
        amp = pd.to_numeric(sub.get("sigmoid_amplitude_mag", pd.Series(dtype=float)), errors="coerce").dropna()
        mean_amp = float(amp.mean()) if not amp.empty else np.nan
        mean_err = _mean_positive(sub.get("extinction_error_mag", pd.Series(dtype=float)))
        mean_sig = _mean_positive(sub.get("stellar_sample_sigma_mag", pd.Series(dtype=float)))
        ratio = mean_amp / mean_err if (np.isfinite(mean_err) and mean_err > 0) else np.nan
        rows.append(
            {
                "id": bid,
                "sb_label": f"SB{bid}",
                "shape": b["shape"],
                "mark": int(b["mark"]),
                "source_type": b["source_type"],
                "fit_std_pc": b["fit_std_pc"],
                "n_shell_components": int(len(sub)),
                "mean_sigmoid_amplitude_mag": mean_amp,
                "amplitude_median_mag": float(amp.median()) if not amp.empty else np.nan,
                "mean_extinction_error_mag": mean_err,
                "three_sigma_error_mag": 3.0 * mean_err if np.isfinite(mean_err) else np.nan,
                "mean_stellar_sample_sigma_mag": mean_sig,
                "amplitude_to_error_ratio": ratio,
                "above_threshold": bool(ratio > AMP_ERROR_RATIO_THRESHOLD) if np.isfinite(ratio) else False,
            }
        )
    return pd.DataFrame(rows)


def summarize_shell_pixel_3sigma(bubbles, detail):
    """Summarize the fraction of unique shell pixels detected at >=3 sigma.

    A shell pixel is one unique (OSB ID, HEALPix pixel) pair containing at
    least one fitted dust component whose center falls inside the adopted
    shell mask. Only component hits with finite, positive extinction errors
    enter the denominator. A pixel passes when at least one of its valid
    shell-associated components satisfies

        sigmoid_amplitude_mag >= AMP_ERROR_RATIO_THRESHOLD * extinction_error_mag.
    """
    base = bubbles[["id"]].copy()
    base["id"] = base["id"].astype(int)

    if detail.empty:
        base["valid_unique_shell_pixels"] = 0
        base["unique_shell_pixels_ge_3sigma"] = 0
        base["unique_shell_pixel_fraction_ge_3sigma"] = np.nan
        return base, {
            "valid_shell_component_hits": 0,
            "invalid_error_component_hits": 0,
            "valid_unique_shell_pixels": 0,
            "unique_shell_pixels_ge_3sigma": 0,
            "unique_shell_pixel_fraction_ge_3sigma": np.nan,
        }

    amplitude = pd.to_numeric(detail["sigmoid_amplitude_mag"], errors="coerce").to_numpy(float)
    error = pd.to_numeric(detail["extinction_error_mag"], errors="coerce").to_numpy(float)
    valid = np.isfinite(amplitude) & np.isfinite(error) & (error > 0.0)

    component_hits = detail.loc[valid, ["id", "healpix_pixel"]].copy()
    component_hits["id"] = component_hits["id"].astype(int)
    component_hits["at_least_3sigma"] = (
        amplitude[valid] >= AMP_ERROR_RATIO_THRESHOLD * error[valid]
    )

    per_pixel = (
        component_hits.groupby(["id", "healpix_pixel"], as_index=False, sort=True)
        .agg(at_least_3sigma=("at_least_3sigma", "any"))
    )
    per_osb = (
        per_pixel.groupby("id", as_index=False, sort=True)
        .agg(
            valid_unique_shell_pixels=("healpix_pixel", "size"),
            unique_shell_pixels_ge_3sigma=("at_least_3sigma", "sum"),
        )
    )
    per_osb["unique_shell_pixel_fraction_ge_3sigma"] = (
        per_osb["unique_shell_pixels_ge_3sigma"] / per_osb["valid_unique_shell_pixels"]
    )

    result = base.merge(per_osb, on="id", how="left")
    for column in ["valid_unique_shell_pixels", "unique_shell_pixels_ge_3sigma"]:
        result[column] = result[column].fillna(0).astype(int)
    result["unique_shell_pixel_fraction_ge_3sigma"] = np.where(
        result["valid_unique_shell_pixels"] > 0,
        result["unique_shell_pixels_ge_3sigma"] / result["valid_unique_shell_pixels"],
        np.nan,
    )

    total_pixels = int(len(per_pixel))
    total_ge_3sigma = int(per_pixel["at_least_3sigma"].sum())
    totals = {
        "valid_shell_component_hits": int(valid.sum()),
        "invalid_error_component_hits": int((~valid).sum()),
        "valid_unique_shell_pixels": total_pixels,
        "unique_shell_pixels_ge_3sigma": total_ge_3sigma,
        "unique_shell_pixel_fraction_ge_3sigma": (
            float(total_ge_3sigma / total_pixels) if total_pixels else np.nan
        ),
    }
    return result, totals


# ----------------------------------------------------------------------------
# Figure styling and layout configuration for reproducible rendering.
# ----------------------------------------------------------------------------
def configure_matplotlib():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.unicode_minus": False,
            "figure.dpi": 140,
            "savefig.dpi": 300,
            "font.size": 17,
            "axes.titlesize": 18,
            "axes.labelsize": 19,
            "axes.linewidth": 1.2,
            "xtick.labelsize": 16,
            "ytick.labelsize": 16,
            "legend.fontsize": 17,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 5,
            "ytick.major.size": 5,
            "xtick.major.width": 1.1,
            "ytick.major.width": 1.1,
        }
    )


def draw(summary, output):
    df = summary.sort_values("id").reset_index(drop=True)
    x = np.arange(len(df))
    amp = df["mean_sigmoid_amplitude_mag"].to_numpy(float)
    err = df["mean_extinction_error_mag"].to_numpy(float)
    three_err = df["three_sigma_error_mag"].to_numpy(float)
    sig = df["mean_stellar_sample_sigma_mag"].to_numpy(float)

    fig, ax = plt.subplots(figsize=(15.8, 8.0))
    w = 0.27
    ax.bar(x - w, amp, width=w, color="#4C78A8", alpha=0.92, label="Mean sigmoid amplitude")
    ax.bar(x, three_err, width=w, color="#D1495B", alpha=0.88, label="3 x mean extinction error")
    ax.bar(x + w, sig, width=w, color="#8E6C8A", alpha=0.88, label="Mean stellar sample sigma")

    median_lines = [
        ("Median amplitude", amp, "#4C78A8"),
        ("Median error", err, "#D1495B"),
        ("Median stellar sigma", sig, "#8E6C8A"),
    ]
    median_labels = []
    for label, values, color in median_lines:
        vals = values[np.isfinite(values)]
        if vals.size:
            median = float(np.median(vals))
            ax.axhline(median, color=color, linestyle="--", linewidth=1.6, alpha=0.75)
            median_labels.append(f"{label}: {median:.4f} mag")

    ax.set_ylabel("Magnitude (mag)")
    ax.set_xlabel("SB ID", labelpad=10)
    ax.set_xticks(x)
    ax.set_xticklabels(df["id"].astype(int).astype(str), rotation=0, ha="center")
    ax.tick_params(axis="both", pad=5)
    ax.grid(axis="y", alpha=0.22, linewidth=0.9)
    ax.legend(frameon=False, ncol=3, loc="upper left", bbox_to_anchor=(0.0, 1.02), columnspacing=1.4)
    if median_labels:
        ax.text(
            0.99,
            0.90,
            "\n".join(median_labels),
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=16,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 5},
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
def build_parser():
    """Build the command-line interface for Script 9."""
    p = argparse.ArgumentParser(description='Run Script 9: OSBs shell amplitude test.')
    p.add_argument("--data-fits", type=Path, default=None, help=argparse.SUPPRESS)
    p.add_argument("--chunk-size", type=int, default=CHUNK_SIZE, help='Command-line option for the documented workflow.')
    p.add_argument("--allow-beyond-max-distance", action="store_true", help='Command-line option for the documented workflow.')
    p.add_argument("--write-detail", action="store_true", help='Input file used by this stage.')
    p.add_argument("--quiet", action="store_true")
    return p


def main():
    """Run Script 9 from validated inputs to the documented outputs."""
    args = build_parser().parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)

    data_fits = resolve_data_fits(args.data_fits)
    bubbles = load_bubbles(GEOM_TABLE)
    n_ell = int((bubbles["shape"] == "ellipsoid").sum())
    n_cyl = int((bubbles["shape"] == "cylinder").sum())
    print(f"Loaded OSB geometry table: {GEOM_TABLE}")
    print(f"Objects to scan: {len(bubbles)} ({n_ell} ellipsoids, {n_cyl} cylinders)")
    print(f"      data_v3.fits = {data_fits}")
    print(f"      extinction-error cube = {ERROR_FITS}")

    print("[2/4] scanning shell components in the dustmaps3d cube ...")
    detail = scan_shell_components(bubbles, data_fits, args.chunk_size, args.allow_beyond_max_distance, args.quiet)

    print("[3/4] attaching extinction errors and stellar-sample sigma ...")
    detail = attach_error_and_sigma(detail, data_fits, ERROR_FITS)

    summary = summarize(bubbles, detail)
    pixel_3sigma, pixel_3sigma_totals = summarize_shell_pixel_3sigma(bubbles, detail)
    summary.to_csv(OUT_SUMMARY, index=False, encoding="utf-8-sig")
    if args.write_detail:
        detail.to_csv(OUT_DETAIL, index=False, encoding="utf-8-sig")

    print("[4/4] plotting ...")
    configure_matplotlib()
    draw(summary, OUT_FIG)

    print(f"summary: {OUT_SUMMARY}")
    if args.write_detail:
        print(f"detail: {OUT_DETAIL}")
    print(f"figure: {OUT_FIG}")
    show = summary[
        ["sb_label", "shape", "n_shell_components", "mean_sigmoid_amplitude_mag",
         "mean_extinction_error_mag", "three_sigma_error_mag",
         "mean_stellar_sample_sigma_mag", "amplitude_to_error_ratio", "above_threshold"]
    ]
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(show.to_string(index=False))
    n_above = int(summary["above_threshold"].sum())
    print(f"OSBs above the shell-amplitude threshold: {n_above}/{len(summary)}")

    pixel_show = pixel_3sigma.copy()
    pixel_show["fraction_ge_3sigma_percent"] = (
        100.0 * pixel_show["unique_shell_pixel_fraction_ge_3sigma"]
    )
    print(
        "\nPixel-level >=3-sigma shell fractions "
        "(one unique OSB-HEALPix pair; any valid shell component may qualify):"
    )
    with pd.option_context("display.max_rows", None, "display.width", 200):
        print(
            pixel_show[
                [
                    "id",
                    "valid_unique_shell_pixels",
                    "unique_shell_pixels_ge_3sigma",
                    "fraction_ge_3sigma_percent",
                ]
            ].to_string(
                index=False,
                formatters={"fraction_ge_3sigma_percent": lambda value: f"{value:.4f}"},
            )
        )
    total_fraction = pixel_3sigma_totals["unique_shell_pixel_fraction_ge_3sigma"]
    total_percent = 100.0 * total_fraction if np.isfinite(total_fraction) else np.nan
    print(
        "Pooled pixel-level >=3-sigma fraction: "
        f"{pixel_3sigma_totals['unique_shell_pixels_ge_3sigma']}/"
        f"{pixel_3sigma_totals['valid_unique_shell_pixels']} "
        f"({total_percent:.4f}%)"
    )
    print(
        "Excluded shell-component hits with non-finite or non-positive "
        f"extinction errors: {pixel_3sigma_totals['invalid_error_component_hits']}"
    )

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
