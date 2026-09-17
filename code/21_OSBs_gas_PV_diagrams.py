"""Script 21: OSBs gas PV diagrams

Purpose
-------
Draw CO/HI longitude-velocity comparisons for three shell velocity prescriptions.

Method overview
---------------
1. Load the compact CO/HI longitude-velocity backgrounds, molecular clouds and
   theoretical/observed spiral-arm tracks.
2. Project each fitted shell into longitude and line-of-sight velocity using the three
   configured kinematic prescriptions.
3. Render the reference panel and model variants with the same survey backgrounds so the
   shell-velocity assumptions can be compared directly.

Main inputs
-----------
- ../data/COdata/MWISP_12CO_clustered_MCs_PV_projection.fits
- ../data/COdata/CfA_12CO_PV_cache.npz
- ../data/HIdata/HI4PI_HI_PV_cache.npy
- ../data/HIdata/HI4PI_HI_PV_axes.npz
- ../data/Reid2019_spirals_xy.csv
- ../data/Cantilever/*_lbvRBD.2019
- ../data/MCs.csv
- ../results/superbubble_final_fit_parameters.csv

Main outputs
------------
- ../results/figures/21_PV_*.png; CO/HI PV comparison figures.

Figure/table role
-----------------
../results/figures/21_PV_*.png; CO/HI PV comparison figures.

Runtime and data notes
----------------------
Reference runtime: 27.98 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.sans-serif"] = [
    "Microsoft YaHei",
    "SimHei",
    "Noto Sans CJK SC",
    "Arial Unicode MS",
    "DejaVu Sans",
]
matplotlib.rcParams["axes.unicode_minus"] = False
matplotlib.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 15,
    "axes.labelsize": 14,
    "xtick.labelsize": 13,
    "ytick.labelsize": 13,
    "legend.fontsize": 14,
    "figure.titlesize": 16,
})

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib import patheffects
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "21_co_hi_PV_diagrams"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"
CO_CACHE_DIR = DATA_DIR / "COdata"
HI_CACHE_DIR = DATA_DIR / "HIdata"
MWISP_OD_FITS = DATA_DIR / "COdata" / "MWISP_12CO_clustered_MCs_PV_projection.fits"
CFA_CACHE = CO_CACHE_DIR / "CfA_12CO_PV_cache.npz"
HI_LV_CACHE = HI_CACHE_DIR / "HI4PI_HI_PV_cache.npy"
HI_AXES = HI_CACHE_DIR / "HI4PI_HI_PV_axes.npz"
FINAL_SB_TABLE = FINAL_DIR / "superbubble_final_fit_parameters.csv"
MCS_CATALOG = DATA_DIR / "MCs.csv"
SPIRAL_XY_PATH = DATA_DIR / "Reid2019_spirals_xy.csv"
CANTILEVER_DIR = DATA_DIR / "Cantilever"

FIG_DIR = FINAL_FIG_DIR
SUMMARY_DIR = OUT_DIR
REFERENCE_PANEL_PNG = FIG_DIR / "21_PV_reid2019_mcs_reference.png"

CENTER_DEG = 180.0
EXPANSION_KMS = 6.0
N_ELLIPSE_SAMPLES = 361
Y_LIMIT = (-50.0, 50.0)
REFERENCE_Y_LIMIT = (-60.0, 60.0)
CO_CMAP = "Purples"
BG_CMAP = "Purples"
HI_BG_CMAP = "Purples"
MAX_BACKGROUND_LON_BINS = 7000
HI_LINEAR_UPPER_PERCENTILE = 99.8
MCS_CLOUD_COLOR = "#ff8c00"
PANEL_TITLE_FONTSIZE = 15
MAIN_TITLE_FONTSIZE = 16
LEGEND_FONTSIZE = 14
OSB_LABEL_FONTSIZE = 6.5
COLORBAR_LABEL_FONTSIZE = 14
COLORBAR_TICK_FONTSIZE = 13
REID_ARM_LINE_EFFECTS = [
    patheffects.Stroke(linewidth=3.0, foreground="white", alpha=0.94),
    patheffects.Normal(),
]
REID_OBS_ARM_LINE_EFFECTS = [
    patheffects.Stroke(linewidth=2.5, foreground="white", alpha=0.9),
    patheffects.Normal(),
]
OVERLAY_COLORS = [
    "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00",
    "#a65628", "#f781bf", "#1b1b1b", "#00a087", "#3c5488",
    "#d62728", "#008fd5", "#2ca02c", "#9467bd", "#ff6f00",
    "#8c564b", "#c2185b", "#005f73", "#6a994e", "#7b2cbf",
    "#bc6c25", "#0077b6", "#d00000", "#2d6a4f", "#7209b7",
    "#f77f00", "#006d77", "#9d0208", "#3a0ca3", "#588157",
]
OVERLAY_LINE_EFFECTS = [
    patheffects.Stroke(linewidth=1.65, foreground="white", alpha=0.88),
    patheffects.Stroke(linewidth=1.05, foreground="black", alpha=0.48),
    patheffects.Normal(),
]

R0 = 8.15
THETA0 = 236.0
U_SUN = 10.7
V_SUN = 11.9
W_SUN = 7.7
A2 = 0.96
A3 = 1.62
U_STD = 10.3
V_STD = 15.3
W_STD = 7.7

FOCUS_THEORY_ARMS = [
    ("Perseus Arm", ["Perseus"], "#000000"),
    ("Local Arm", ["Local"], "#17becf"),
    ("Sagittarius-Carina Arm", ["Sagittarius", "SgF", "Carina", "CrF"], "#a03fd1"),
    ("Local Arm Spur", ["LoS"], "#2ca02c"),
]

FOCUS_OBS_ARMS = [
    ("Perseus Arm", ["Per"], "#000000"),
    ("Local Arm", ["Loc", "AqS", "AqR"], "#17becf"),
    ("Sagittarius-Carina Arm", ["SgN", "SgF", "CrN", "CrF"], "#a03fd1"),
    ("Local Arm Spur", ["LoS"], "#2ca02c"),
]

VARIANTS = {
    "center": {
        "out_png": FIG_DIR / "21_PV_model1_center_velocity_plus6.png",
        "out_csv": SUMMARY_DIR / "21_PV_model1_center_velocity_plus6_summary.csv",
        "label": "center circular velocity + 6 km/s radial expansion",
        "title": "Model 1: shell inherits superbubble-center circular velocity, then expands at 6 km/s",
        "summary_tag": "center_velocity_plus6",
    },
    "local": {
        "out_png": FIG_DIR / "21_PV_model2_local_circular_plus6.png",
        "out_csv": SUMMARY_DIR / "21_PV_model2_local_circular_plus6_summary.csv",
        "label": "local shell circular velocity + 6 km/s radial expansion",
        "title": "Model 2: shell is fully circularized by Galactic rotation, then expands at 6 km/s",
        "summary_tag": "local_circular_plus6",
    },
    "halfscale": {
        "out_png": FIG_DIR / "21_PV_model3_halfscale_velocity_plus6.png",
        "out_csv": SUMMARY_DIR / "21_PV_model3_halfscale_velocity_plus6_summary.csv",
        "label": "original positions, half-scale circular velocity reference + 6 km/s radial expansion",
        "title": "Model 3: swept-gas mixed velocity from half-scale xy ellipse, then expands at 6 km/s",
        "summary_tag": "halfscale_velocity_plus6",
    },
}


def reid2019_urc(r_gal: np.ndarray) -> np.ndarray:
    r = np.asarray(r_gal, dtype=float)
    lam = (A3 / 1.5) ** 5
    r_opt = A2 * R0
    rho = np.clip(r / r_opt, 1e-6, None)
    log_lam = np.log10(lam)
    term1 = 200.0 * lam**0.41
    top = 0.75 * np.exp(-0.4 * lam)
    bot = 0.47 + 2.25 * lam**0.4
    term2 = np.sqrt(0.80 + 0.49 * log_lam + (top / bot))
    top = 1.97 * rho**1.22
    bot = (rho**2 + 0.61) ** 1.43
    term3 = (0.72 + 0.44 * log_lam) * (top / bot)
    top = rho**2
    bot = rho**2 + 2.25 * lam**0.4
    term4 = 1.60 * np.exp(-0.4 * lam) * (top / bot)
    return (term1 / term2) * np.sqrt(term3 + term4)


def lbd_to_vlsr(l_deg: np.ndarray, b_deg: np.ndarray, d: np.ndarray) -> np.ndarray:
    l_rad = np.radians(np.asarray(l_deg, dtype=float))
    b_rad = np.radians(np.asarray(b_deg, dtype=float))
    d = np.asarray(d, dtype=float)
    d_proj = d * np.cos(b_rad)
    x_gc = d_proj * np.sin(l_rad)
    y_gc = R0 - d_proj * np.cos(l_rad)
    r_gal = np.sqrt(x_gc**2 + y_gc**2)
    sin_b = np.sin(b_rad)
    cos_b = np.cos(b_rad)
    valid = (r_gal > 0.2) & np.isfinite(d) & (d > 0.0)
    theta_r = reid2019_urc(r_gal)
    r_hat_x = np.divide(x_gc, r_gal, out=np.zeros_like(x_gc), where=r_gal > 0.0)
    r_hat_y = np.divide(y_gc, r_gal, out=np.zeros_like(y_gc), where=r_gal > 0.0)
    t_hat_x = r_hat_y
    t_hat_y = -r_hat_x
    v_src_x = theta_r * t_hat_x
    v_src_y = theta_r * t_hat_y
    v_src_z = np.zeros_like(v_src_x)
    v_sun_x = THETA0 + V_SUN
    v_sun_y = -U_SUN
    v_sun_z = W_SUN
    los_x = np.cos(b_rad) * np.sin(l_rad)
    los_y = -np.cos(b_rad) * np.cos(l_rad)
    los_z = np.sin(b_rad)
    v_helio = ((v_src_x - v_sun_x) * los_x + (v_src_y - v_sun_y) * los_y + (v_src_z - v_sun_z) * los_z)
    v_lsr = v_helio + (U_STD * np.cos(l_rad) + V_STD * np.sin(l_rad)) * cos_b + W_STD * sin_b
    v_lsr[~valid] = np.nan
    return v_lsr


def plot_longitude_for_center(l_deg: np.ndarray, center_deg: float = 0.0) -> np.ndarray:
    return np.mod(np.asarray(l_deg, dtype=float) - center_deg + 180.0, 360.0)


def absolute_longitude_ticklabels(center_deg: float):
    labels = []
    for tick in np.arange(0, 361, 60):
        val = int(round((tick + center_deg - 180.0) % 360.0))
        if tick == 360 and val == 0:
            val = 360
        labels.append(str(val))
    return labels


def xy_to_lv(x: np.ndarray, y: np.ndarray):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    dx = x
    dy = y - R0
    d = np.sqrt(dx**2 + dy**2)
    l_deg = (np.degrees(np.arctan2(dx, -dy)) + 360.0) % 360.0
    b_deg = np.zeros_like(l_deg)
    v_lsr = lbd_to_vlsr(l_deg, b_deg, d)
    return l_deg, v_lsr


def split_original_path(x: np.ndarray, y: np.ndarray, l_deg: np.ndarray, v_lsr: np.ndarray):
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(l_deg) & np.isfinite(v_lsr)
    if not np.any(valid):
        return []
    x = np.asarray(x[valid], dtype=float)
    y = np.asarray(y[valid], dtype=float)
    l_deg = np.asarray(l_deg[valid], dtype=float)
    v_lsr = np.asarray(v_lsr[valid], dtype=float)
    ds = np.hypot(np.diff(x), np.diff(y))
    if ds.size == 0:
        return [(l_deg, v_lsr)]
    positive = ds[ds > 0]
    gap_thresh = 0.5 if positive.size == 0 else max(0.35, float(np.median(positive) * 3.5), float(np.percentile(positive, 95)))
    split_idx = np.where((ds > gap_thresh) | (np.abs(np.diff(l_deg)) > 180.0))[0] + 1
    return list(zip(np.split(l_deg, split_idx), np.split(v_lsr, split_idx)))


def split_centered_segments(l_deg: np.ndarray, v_lsr: np.ndarray, center_deg: float = 0.0):
    valid = np.isfinite(l_deg) & np.isfinite(v_lsr)
    if not np.any(valid):
        return []
    l_plot = plot_longitude_for_center(np.asarray(l_deg[valid], dtype=float), center_deg=center_deg)
    v_plot = np.asarray(v_lsr[valid], dtype=float)
    if l_plot.size < 2:
        return [(l_plot, v_plot)]
    split_idx = np.where(np.abs(np.diff(l_plot)) > 300.0)[0] + 1
    return list(zip(np.split(l_plot, split_idx), np.split(v_plot, split_idx)))


def read_cantilever_tracks():
    code_to_data = {}
    for file in sorted(CANTILEVER_DIR.glob("*_lbvRBD.2019")):
        code = file.stem.split("_")[0]
        df = pd.read_csv(file, sep=r"\s+", comment="!", header=None, names=["l", "b", "v", "R", "B", "D"])
        if not df.empty:
            code_to_data[code] = df
    return code_to_data


def gc_velocity_to_vlsr(base, l_deg: np.ndarray, b_deg: np.ndarray, vx_gc, vy_gc, vz_gc) -> np.ndarray:
    l_rad = np.radians(l_deg)
    b_rad = np.radians(b_deg)
    sin_b = np.sin(b_rad)
    cos_b = np.cos(b_rad)
    los_x = cos_b * np.sin(l_rad)
    los_y = -cos_b * np.cos(l_rad)
    los_z = sin_b
    v_sun_x = base.THETA0 + base.V_SUN
    v_sun_y = -base.U_SUN
    v_sun_z = base.W_SUN
    v_helio = (
        (vx_gc - v_sun_x) * los_x
        + (vy_gc - v_sun_y) * los_y
        + (vz_gc - v_sun_z) * los_z
    )
    return (
        v_helio
        + (base.U_STD * np.cos(l_rad) + base.V_STD * np.sin(l_rad)) * cos_b
        + base.W_STD * sin_b
    )


def circular_velocity_vector_at_heliocentric_xy(base, x_helio: float, y_helio: float) -> tuple[float, float, float]:
    x_gc = float(y_helio)
    y_gc = base.R0 - float(x_helio)
    r_gal = np.hypot(x_gc, y_gc)
    if not np.isfinite(r_gal) or r_gal <= 0.2:
        return np.nan, np.nan, np.nan
    theta_r = float(base.reid2019_urc(np.array([r_gal]))[0])
    r_hat_x = x_gc / r_gal
    r_hat_y = y_gc / r_gal
    return theta_r * r_hat_y, -theta_r * r_hat_x, 0.0


class LocalBase:
    CENTER_DEG = CENTER_DEG
    R0 = R0
    THETA0 = THETA0
    U_SUN = U_SUN
    V_SUN = V_SUN
    W_SUN = W_SUN
    U_STD = U_STD
    V_STD = V_STD
    W_STD = W_STD
    SPIRAL_XY_PATH = SPIRAL_XY_PATH
    FOCUS_THEORY_ARMS = FOCUS_THEORY_ARMS
    FOCUS_OBS_ARMS = FOCUS_OBS_ARMS
    reid2019_urc = staticmethod(reid2019_urc)
    lbd_to_vlsr = staticmethod(lbd_to_vlsr)
    plot_longitude_for_center = staticmethod(plot_longitude_for_center)
    absolute_longitude_ticklabels = staticmethod(absolute_longitude_ticklabels)
    xy_to_lv = staticmethod(xy_to_lv)
    split_original_path = staticmethod(split_original_path)
    split_centered_segments = staticmethod(split_centered_segments)
    read_cantilever_tracks = staticmethod(read_cantilever_tracks)


class LocalModel:
    gc_velocity_to_vlsr = staticmethod(gc_velocity_to_vlsr)
    circular_velocity_vector_at_heliocentric_xy = staticmethod(circular_velocity_vector_at_heliocentric_xy)


def axis_edges(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    step = float(np.median(np.diff(values))) if values.size > 1 else 1.0
    edges = np.empty(values.size + 1, dtype=np.float64)
    edges[1:-1] = 0.5 * (values[:-1] + values[1:])
    edges[0] = values[0] - 0.5 * step
    edges[-1] = values[-1] + 0.5 * step
    return edges


def robust_log_limits(data: np.ndarray) -> tuple[float, float]:
    positive = np.asarray(data)[np.isfinite(data) & (data > 0.0)]
    if positive.size == 0:
        return 1e-3, 1.0
    vmin = max(float(np.percentile(positive, 1.5)), 1e-4)
    vmax = float(np.percentile(positive, 99.8))
    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = max(vmin * 10.0, float(np.nanmax(positive)))
    return vmin, vmax


def percentile_limits(image: np.ndarray, low: float = 1.0, high: float = 99.5) -> tuple[float, float]:
    finite = np.asarray(image)[np.isfinite(image)]
    if finite.size == 0:
        return 0.0, 1.0
    vmin, vmax = np.percentile(finite, [low, high])
    if np.isclose(vmin, vmax):
        vmax = vmin + 1.0
    return float(vmin), float(vmax)


def robust_linear_limits(data: np.ndarray, upper_pct: float = 95.0) -> tuple[float, float]:
    finite = np.asarray(data)[np.isfinite(data)]
    positive = finite[finite > 0.0]
    if positive.size == 0:
        return 0.0, 1.0
    vmax = float(np.percentile(positive, upper_pct))
    if not np.isfinite(vmax) or vmax <= 0.0:
        vmax = float(np.nanmax(positive))
    return 0.0, max(vmax, 1.0)


def downsample_longitude(
    lon_plot: np.ndarray,
    lv_map: np.ndarray,
    max_bins: int = MAX_BACKGROUND_LON_BINS,
) -> tuple[np.ndarray, np.ndarray]:
    if lv_map.shape[1] <= max_bins:
        return lon_plot, lv_map

    factor = int(np.ceil(lv_map.shape[1] / max_bins))
    n = (lv_map.shape[1] // factor) * factor
    lon_trim = lon_plot[:n].reshape(-1, factor)
    lv_trim = lv_map[:, :n].reshape(lv_map.shape[0], -1, factor)
    lon_ds = np.nanmean(lon_trim, axis=1)
    lv_ds = np.nanmean(lv_trim, axis=2)

    if n < lv_map.shape[1]:
        lon_tail = np.array([np.nanmean(lon_plot[n:])])
        lv_tail = np.nanmean(lv_map[:, n:], axis=1, keepdims=True)
        lon_ds = np.concatenate([lon_ds, lon_tail])
        lv_ds = np.concatenate([lv_ds, lv_tail], axis=1)
    return lon_ds, lv_ds


def load_mwisp_od_fits_for_plot(base):
    from astropy.io import fits

    with fits.open(MWISP_OD_FITS, memmap=True) as hdul:
        header = hdul[0].header
        data = np.asarray(hdul[0].data, dtype=np.float32)

    nlon = data.shape[1]
    nvel = data.shape[0]
    lon = (
        (np.arange(nlon, dtype=np.float64) + 1.0 - float(header["CRPIX1"]))
        * float(header["CDELT1"])
        + float(header["CRVAL1"])
    )
    lon = np.mod(lon, 360.0)
    vel = (
        (np.arange(nvel, dtype=np.float64) + 1.0 - float(header.get("CRPIX3", 1.0)))
        * float(header.get("CDELT3", 1000.0))
        + float(header.get("CRVAL3", 0.0))
    ) / 1000.0

    vel_mask = (vel >= Y_LIMIT[0] - 8.0) & (vel <= Y_LIMIT[1] + 8.0)
    vel = vel[vel_mask]
    lv_map = data[vel_mask, :]

    lon_plot = base.plot_longitude_for_center(lon, center_deg=CENTER_DEG)
    order = np.argsort(lon_plot)
    lon_sorted = lon_plot[order]
    lv_sorted = lv_map[:, order]
    lon_sorted, lv_sorted = downsample_longitude(lon_sorted, lv_sorted)

    disp = np.ma.array(np.array(lv_sorted, copy=True), mask=(~np.isfinite(lv_sorted)) | (lv_sorted <= 0.0))
    vmin, vmax = robust_log_limits(lv_sorted)
    return axis_edges(lon_sorted), axis_edges(vel), disp, vmin, vmax, "OD_12CO masked LV"


def load_lv_cache_for_plot(cache_path: Path, base):
    with np.load(cache_path) as cached:
        lv_map = cached["lv_map"].astype(np.float32, copy=False)
        lon = cached["lon"].astype(np.float64, copy=False)
        vel = cached["vel"].astype(np.float64, copy=False)
        if "b_limit" in cached.files:
            b_text = f"|b| <= {float(cached['b_limit'][0]):.0f} deg"
        else:
            b_text = f"b={float(cached['bmin'][0]):.2f} to {float(cached['bmax'][0]):.2f} deg"

    lon_plot = base.plot_longitude_for_center(lon, center_deg=CENTER_DEG)
    order = np.argsort(lon_plot)
    lon_sorted = lon_plot[order]
    lv_sorted = lv_map[:, order]
    disp = np.ma.array(np.array(lv_sorted, copy=True), mask=(~np.isfinite(lv_sorted)) | (lv_sorted <= 0.0))
    vmin, vmax = robust_log_limits(lv_sorted)
    return axis_edges(lon_sorted), axis_edges(vel), disp, vmin, vmax, b_text


def draw_lv_background(ax, cache_path: Path, base, alpha: float = 0.86):
    lon_edges, vel_edges, disp, vmin, vmax, b_text = load_lv_cache_for_plot(cache_path, base)
    cmap = plt.get_cmap(BG_CMAP).copy()
    cmap.set_bad((1, 1, 1, 0))
    im = ax.pcolormesh(
        lon_edges,
        vel_edges,
        disp,
        cmap=cmap,
        norm=LogNorm(vmin=vmin, vmax=vmax),
        shading="auto",
        alpha=alpha,
        zorder=0,
    )
    return im, b_text, vmin, vmax


def draw_mwisp_od_background(ax, base, alpha: float = 0.86):
    lon_edges, vel_edges, disp, vmin, vmax, b_text = load_mwisp_od_fits_for_plot(base)
    cmap = plt.get_cmap(BG_CMAP).copy()
    cmap.set_bad((1, 1, 1, 0))
    im = ax.pcolormesh(
        lon_edges,
        vel_edges,
        disp,
        cmap=cmap,
        norm=LogNorm(vmin=vmin, vmax=vmax),
        shading="auto",
        alpha=alpha,
        zorder=0,
    )
    return im, b_text, vmin, vmax


def load_hi_lv_cache_for_plot(base):
    lv_map = np.load(HI_LV_CACHE).astype(np.float32, copy=False)
    with np.load(HI_AXES) as axes:
        lon = axes["lon_deg"].astype(np.float64, copy=False)
        vel = axes["vel_kms"].astype(np.float64, copy=False)

    lon_plot = base.plot_longitude_for_center(lon, center_deg=CENTER_DEG)
    order = np.argsort(lon_plot)
    lon_sorted = lon_plot[order]
    lv_sorted = lv_map[:, order]
    disp = np.ma.array(np.array(lv_sorted, copy=True), mask=~np.isfinite(lv_sorted))

    vmin, vmax = percentile_limits(lv_sorted, low=3.0, high=99.7)
    return axis_edges(lon_sorted), axis_edges(vel), disp, vmin, vmax, "|b| <= 15 deg"


def draw_hi_background(ax, base, alpha: float = 1.00):
    lon_edges, vel_edges, disp, vmin, vmax, b_text = load_hi_lv_cache_for_plot(base)
    cmap = plt.get_cmap(HI_BG_CMAP).copy()
    cmap.set_bad("#f4f4f4")
    im = ax.pcolormesh(
        lon_edges,
        vel_edges,
        disp,
        cmap=cmap,
        norm=Normalize(vmin=vmin, vmax=vmax, clip=True),
        shading="auto",
        alpha=alpha,
        zorder=0,
    )
    return im, b_text, vmin, vmax


def plot_reid2019_theory_arms(ax, base) -> list[Line2D]:
    spiral_xy = pd.read_csv(base.SPIRAL_XY_PATH)
    handles = []
    for name, codes, color in base.FOCUS_THEORY_ARMS:
        handle = Line2D([0], [0], color=color, lw=2.1, label=name)
        handles.append(handle)
        for code in codes:
            arm_df = spiral_xy[spiral_xy["arm"] == code]
            if arm_df.empty:
                continue
            x = arm_df["xx"].to_numpy(dtype=float)
            y = arm_df["yy"].to_numpy(dtype=float)
            l_deg, v_lsr = base.xy_to_lv(x, y)
            for ll, vv in base.split_original_path(x, y, l_deg, v_lsr):
                for llc, vvc in base.split_centered_segments(ll, vv, center_deg=CENTER_DEG):
                    if llc.size < 2:
                        continue
                    line, = ax.plot(
                        llc,
                        vvc,
                        color=color,
                        linewidth=1.55,
                        alpha=0.97,
                        zorder=7,
                        path_effects=REID_ARM_LINE_EFFECTS,
                    )
    return handles


def plot_reid2019_observed_arms(ax, base) -> int:
    cantilever = base.read_cantilever_tracks()
    n_tracks = 0
    for _, codes, color in base.FOCUS_OBS_ARMS:
        for code in codes:
            if code not in cantilever:
                continue
            df = cantilever[code]
            l_deg = df["l"].to_numpy(dtype=float)
            v_lsr = df["v"].to_numpy(dtype=float)
            for llc, vvc in base.split_centered_segments(l_deg, v_lsr, center_deg=CENTER_DEG):
                if llc.size < 2:
                    continue
                line, = ax.plot(
                    llc,
                    vvc,
                    color=color,
                    linewidth=1.2,
                    linestyle="--",
                    dashes=(4, 2),
                    alpha=0.96,
                    zorder=7,
                    path_effects=REID_OBS_ARM_LINE_EFFECTS,
                )
            n_tracks += 1
    return n_tracks


def read_mcs_catalog() -> pd.DataFrame:
    columns = ["Seq", "GLON", "GLAT", "Dist", "Rad", "Area", "Mass", "Sigma", "rho", "Flag"]
    df = pd.read_csv(MCS_CATALOG)
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"molecular-cloud catalogue is missing required columns: {', '.join(missing)}")
    for col in columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def plot_mcs_clouds(ax, base) -> int:
    df = read_mcs_catalog()
    df = df[np.abs(df["GLAT"]) < 5.0].copy()
    l_deg = df["GLON"].to_numpy(dtype=float)
    b_deg = df["GLAT"].to_numpy(dtype=float)
    dist = df["Dist"].to_numpy(dtype=float)
    v_lsr = base.lbd_to_vlsr(l_deg, b_deg, dist)
    l_plot = base.plot_longitude_for_center(l_deg, center_deg=CENTER_DEG)
    mask = np.isfinite(l_plot) & np.isfinite(v_lsr)
    sc = ax.scatter(
        l_plot[mask],
        v_lsr[mask],
        s=14,
        c=MCS_CLOUD_COLOR,
        edgecolors="0.08",
        linewidths=0.22,
        alpha=0.78,
        marker="o",
        zorder=6,
    )
    sc.set_path_effects([
        patheffects.Stroke(linewidth=0.65, foreground="white", alpha=0.76),
        patheffects.Normal(),
    ])
    return int(mask.sum())


def xy_compromise_scale(df: pd.DataFrame, frac: float = 0.5) -> np.ndarray:
    """Scale max-opening xy_plane axes to the mid-cap compromise ellipse."""
    shape = df["shape"].astype(str).str.lower().to_numpy()
    cz = df["center_z_kpc"].to_numpy(float)
    c = df["c_radius_kpc"].to_numpy(float)
    z_base = df["xy_plane_z_kpc"].to_numpy(float)
    mark = df["mark"].to_numpy(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        k_base = np.sqrt(np.clip(1.0 - ((z_base - cz) / c) ** 2, 0.0, None))
        z_apex = np.where(mark == 1, cz - c, cz + c)
        z_mid = z_base + frac * (z_apex - z_base)
        k_mid = np.sqrt(np.clip(1.0 - ((z_mid - cz) / c) ** 2, 0.0, None))
        scale = np.divide(k_mid, k_base, out=np.ones_like(k_mid), where=k_base > 0)
    return np.where((shape == "ellipsoid") & (c > 0), scale, 1.0)


def read_superbubble_table() -> pd.DataFrame:
    needed = [
        "id",
        "mark",
        "shape",
        "xy_model",
        "center_x_kpc",
        "center_y_kpc",
        "center_z_kpc",
        "c_radius_kpc",
        "xy_plane_center_x_kpc",
        "xy_plane_center_y_kpc",
        "xy_plane_z_kpc",
        "xy_plane_a_kpc",
        "xy_plane_b_kpc",
        "xy_plane_angle_deg",
        "xy_plane_definition",
    ]
    source_path = FINAL_SB_TABLE
    if not source_path.exists():
        raise FileNotFoundError(
            f"Could not find the final xy-plane superbubble table: {FINAL_SB_TABLE}"
        )
    df = pd.read_csv(source_path)
    out = df[needed].copy()
    out["xy_plane_a_maxopen_kpc"] = out["xy_plane_a_kpc"].astype(float)
    out["xy_plane_b_maxopen_kpc"] = out["xy_plane_b_kpc"].astype(float)
    out["xy_plane_section_scale"] = xy_compromise_scale(out)
    out["xy_plane_a_compromise_kpc"] = out["xy_plane_a_maxopen_kpc"] * out["xy_plane_section_scale"]
    out["xy_plane_b_compromise_kpc"] = out["xy_plane_b_maxopen_kpc"] * out["xy_plane_section_scale"]
    out["xy_plane_a_kpc"] = out["xy_plane_a_compromise_kpc"]
    out["xy_plane_b_kpc"] = out["xy_plane_b_compromise_kpc"]
    valid = (
        np.isfinite(out["xy_plane_center_x_kpc"])
        & np.isfinite(out["xy_plane_center_y_kpc"])
        & np.isfinite(out["xy_plane_a_kpc"])
        & np.isfinite(out["xy_plane_b_kpc"])
    )
    return out[valid].reset_index(drop=True)


def assign_colors(sb_df: pd.DataFrame) -> dict[int, object]:
    return {
        int(row["id"]): OVERLAY_COLORS[i % len(OVERLAY_COLORS)]
        for i, (_, row) in enumerate(sb_df.iterrows())
    }


def ellipse_points(row: pd.Series, scale: float = 1.0):
    theta = np.linspace(0.0, 2.0 * np.pi, N_ELLIPSE_SAMPLES)
    xc = float(row["xy_plane_center_x_kpc"])
    yc = float(row["xy_plane_center_y_kpc"])
    zc = float(row["xy_plane_z_kpc"])
    a = float(row["xy_plane_a_kpc"]) * float(scale)
    b = float(row["xy_plane_b_kpc"]) * float(scale)
    ang = np.radians(float(row["xy_plane_angle_deg"]))
    x = xc + a * np.cos(theta) * np.cos(ang) - b * np.sin(theta) * np.sin(ang)
    y = yc + a * np.cos(theta) * np.sin(ang) + b * np.sin(theta) * np.cos(ang)
    z = np.full_like(x, zc)
    return x, y, z


def xyz_to_lbd(x: np.ndarray, y: np.ndarray, z: np.ndarray):
    d = np.sqrt(x**2 + y**2 + z**2)
    l_deg = (np.degrees(np.arctan2(y, x)) + 360.0) % 360.0
    b_deg = np.degrees(np.arctan2(z, np.hypot(x, y)))
    invalid = (~np.isfinite(d)) | (d <= 0.0)
    l_deg[invalid] = np.nan
    b_deg[invalid] = np.nan
    d[invalid] = np.nan
    return l_deg, b_deg, d


def circular_vectors_for_points(model, base, x: np.ndarray, y: np.ndarray):
    vx = np.empty_like(x, dtype=float)
    vy = np.empty_like(y, dtype=float)
    vz = np.zeros_like(x, dtype=float)
    for i in range(x.size):
        vx[i], vy[i], vz[i] = model.circular_velocity_vector_at_heliocentric_xy(base, float(x[i]), float(y[i]))
    return vx, vy, vz


def radial_expansion_vector_xy(x: np.ndarray, y: np.ndarray, xc: float, yc: float):
    dx = np.asarray(x, dtype=float) - float(xc)
    dy = np.asarray(y, dtype=float) - float(yc)
    dr = np.hypot(dx, dy)
    ux_h = np.divide(dx, dr, out=np.zeros_like(dx), where=dr > 0.0)
    uy_h = np.divide(dy, dr, out=np.zeros_like(dy), where=dr > 0.0)
    return EXPANSION_KMS * uy_h, -EXPANSION_KMS * ux_h, np.zeros_like(ux_h)


def project_gc_velocity(model, base, l_deg, b_deg, vx, vy, vz):
    v = model.gc_velocity_to_vlsr(base, l_deg, b_deg, vx, vy, vz)
    v[(~np.isfinite(l_deg)) | (~np.isfinite(b_deg))] = np.nan
    return v


def compute_dashed_velocity(variant: str, model, base, row: pd.Series, x, y, z, l_deg, b_deg):
    xc = float(row["xy_plane_center_x_kpc"])
    yc = float(row["xy_plane_center_y_kpc"])

    if variant == "center":
        vx0, vy0, vz0 = model.circular_velocity_vector_at_heliocentric_xy(base, xc, yc)
        vx_base = np.full_like(x, vx0, dtype=float)
        vy_base = np.full_like(y, vy0, dtype=float)
        vz_base = np.full_like(z, vz0, dtype=float)
        velocity_ref_scale = 0.0
    elif variant == "local":
        vx_base, vy_base, vz_base = circular_vectors_for_points(model, base, x, y)
        velocity_ref_scale = 1.0
    elif variant == "halfscale":
        x_ref, y_ref, _ = ellipse_points(row, scale=0.5)
        vx_base, vy_base, vz_base = circular_vectors_for_points(model, base, x_ref, y_ref)
        velocity_ref_scale = 0.5
    else:
        raise ValueError(f"Unknown variant: {variant}")

    vx_exp, vy_exp, vz_exp = radial_expansion_vector_xy(x, y, xc, yc)
    v_dashed = project_gc_velocity(
        model,
        base,
        l_deg,
        b_deg,
        vx_base + vx_exp,
        vy_base + vy_exp,
        vz_base + vz_exp,
    )
    return v_dashed, velocity_ref_scale


def build_overlay_models(variant: str, model, base, sb_df: pd.DataFrame):
    id_to_color = assign_colors(sb_df)
    models = []
    rows = []

    for _, row in sb_df.iterrows():
        sb_id = int(row["id"])
        x, y, z = ellipse_points(row, scale=1.0)
        l_deg, b_deg, d = xyz_to_lbd(x, y, z)
        v_solid = base.lbd_to_vlsr(l_deg, b_deg, d)
        v_dashed, velocity_ref_scale = compute_dashed_velocity(variant, model, base, row, x, y, z, l_deg, b_deg)

        xc = float(row["xy_plane_center_x_kpc"])
        yc = float(row["xy_plane_center_y_kpc"])
        zc = float(row["xy_plane_z_kpc"])
        l_center, b_center, d_center = xyz_to_lbd(np.array([xc]), np.array([yc]), np.array([zc]))
        v_center = float(base.lbd_to_vlsr(l_center, b_center, d_center)[0])
        l_center_plot = float(base.plot_longitude_for_center(l_center, center_deg=CENTER_DEG)[0])

        models.append(
            {
                "id": sb_id,
                "color": id_to_color[sb_id],
                "segments_solid": base.split_original_path(x, y, l_deg, v_solid),
                "segments_dashed": base.split_original_path(x, y, l_deg, v_dashed),
                "l_center_plot": l_center_plot,
                "v_center": v_center,
            }
        )
        rows.append(
            {
                "id": sb_id,
                "mark": int(row["mark"]),
                "shape": str(row["shape"]),
                "xy_model": str(row["xy_model"]),
                "xy_plane_center_x_kpc": xc,
                "xy_plane_center_y_kpc": yc,
                "xy_plane_z_kpc": zc,
                "xy_plane_a_maxopen_kpc": float(row["xy_plane_a_maxopen_kpc"]),
                "xy_plane_b_maxopen_kpc": float(row["xy_plane_b_maxopen_kpc"]),
                "xy_plane_section_scale": float(row["xy_plane_section_scale"]),
                "xy_plane_a_compromise_kpc": float(row["xy_plane_a_compromise_kpc"]),
                "xy_plane_b_compromise_kpc": float(row["xy_plane_b_compromise_kpc"]),
                "xy_plane_a_kpc": float(row["xy_plane_a_kpc"]),
                "xy_plane_b_kpc": float(row["xy_plane_b_kpc"]),
                "xy_plane_angle_deg": float(row["xy_plane_angle_deg"]),
                "xy_plane_definition": str(row["xy_plane_definition"]),
                "solid_velocity_model": "local_shell_circular_velocity",
                "dashed_velocity_model": VARIANTS[variant]["summary_tag"],
                "dashed_position_scale": 1.0,
                "dashed_velocity_reference_scale": velocity_ref_scale,
                "v_expansion_kms": EXPANSION_KMS,
                "l_center_deg": float(l_center[0]),
                "b_center_deg": float(b_center[0]),
                "d_center_kpc": float(d_center[0]),
                "v_center_model_kms": v_center,
                "v_solid_min_kms": float(np.nanmin(v_solid)),
                "v_solid_max_kms": float(np.nanmax(v_solid)),
                "v_dashed_min_kms": float(np.nanmin(v_dashed)),
                "v_dashed_max_kms": float(np.nanmax(v_dashed)),
            }
        )

    return models, pd.DataFrame(rows)


def draw_overlay_models(ax, models):
    for item in models:
        color = item["color"]
        for ll, vv in item["segments_solid"]:
            ax.plot(
                ll,
                vv,
                color=color,
                lw=0.72,
                alpha=1.0,
                linestyle="-",
                zorder=4,
                path_effects=OVERLAY_LINE_EFFECTS,
            )
        for ll, vv in item["segments_dashed"]:
            ax.plot(
                ll,
                vv,
                color=color,
                lw=0.78,
                alpha=1.0,
                linestyle="--",
                dashes=(4, 2),
                zorder=5,
                path_effects=OVERLAY_LINE_EFFECTS,
            )
        if np.isfinite(item["l_center_plot"]) and np.isfinite(item["v_center"]) and Y_LIMIT[0] <= item["v_center"] <= Y_LIMIT[1]:
            ax.text(
                item["l_center_plot"],
                item["v_center"],
                f"SB{item['id']}",
                color=color,
                fontsize=OSB_LABEL_FONTSIZE,
                fontweight="bold",
                ha="center",
                va="center",
                zorder=8,
                path_effects=[patheffects.Stroke(linewidth=1.4, foreground="white", alpha=0.9), patheffects.Normal()],
            )


def style_axis(base, ax, title: str):
    ax.axhline(0.0, color="0.45", lw=0.7, alpha=0.5, zorder=1)
    ax.set_xlim(360.0, 0.0)
    ax.set_ylim(*Y_LIMIT)
    ax.set_xticks(np.arange(0, 361, 60))
    ax.set_xticklabels(base.absolute_longitude_ticklabels(CENTER_DEG))
    ax.grid(color="0.86", linestyle="--", linewidth=0.45, alpha=0.9)
    ax.set_title(title, fontsize=PANEL_TITLE_FONTSIZE)
    ax.set_ylabel(r"$V_{\rm LSR}$ (km s$^{-1}$)")


def render_reference_panel() -> Path:
    REFERENCE_PANEL_PNG.parent.mkdir(parents=True, exist_ok=True)

    base = LocalBase
    base.CENTER_DEG = CENTER_DEG

    fig, ax = plt.subplots(1, 1, figsize=(13.5, 5.8), dpi=360)
    _, mwisp_b_text, mwisp_vmin, mwisp_vmax = draw_mwisp_od_background(ax, base, alpha=0.78)
    n_mcs_clouds = plot_mcs_clouds(ax, base)
    plot_reid2019_theory_arms(ax, base)
    plot_reid2019_observed_arms(ax, base)

    style_axis(
        base,
        ax,
        r"MWISP $^{12}\mathrm{CO}$",
    )
    ax.set_ylim(*REFERENCE_Y_LIMIT)
    ax.set_yticks(np.arange(-60, 61, 20))
    ax.set_xlabel("Galactic Longitude l (deg)", labelpad=8)

    cloud_handle = Line2D(
        [0],
        [0],
        color=MCS_CLOUD_COLOR,
        marker="o",
        linestyle="None",
        markersize=4.5,
        markeredgecolor="0.08",
        label="Molecular clouds (|b| < 5 deg)",
    )
    legend_handles = []
    for arm_name, _, color in FOCUS_THEORY_ARMS:
        if arm_name == "Local Arm":
            continue
        legend_handles.append(
            Line2D([0], [0], color=color, lw=1.6, linestyle="-", label=f"Reid et al. (2019) {arm_name}")
        )
        legend_handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                lw=1.4,
                linestyle="--",
                dashes=(4, 2),
                label=f"Gas-traced {arm_name}",
            )
        )
    legend_handles.append(cloud_handle)
    ax.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.31),
        fontsize=LEGEND_FONTSIZE,
        ncol=2,
        frameon=True,
        framealpha=0.9,
    )

    sm_co = plt.cm.ScalarMappable(norm=LogNorm(vmin=mwisp_vmin, vmax=mwisp_vmax), cmap=plt.get_cmap(BG_CMAP))
    sm_co.set_array([])
    cbar_co = fig.colorbar(sm_co, ax=ax, pad=0.025)
    cbar_co.set_label(r"$\int T_B\,db$ [K deg]", fontsize=COLORBAR_LABEL_FONTSIZE)
    cbar_co.ax.tick_params(labelsize=COLORBAR_TICK_FONTSIZE)

    fig.subplots_adjust(left=0.075, right=0.90, bottom=0.60, top=0.90)
    fig.savefig(REFERENCE_PANEL_PNG, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved model-independent reference panel to: {REFERENCE_PANEL_PNG}")
    return REFERENCE_PANEL_PNG


def render_variant(variant: str = "center") -> Path:
    if variant not in VARIANTS:
        raise ValueError(f"variant must be one of {sorted(VARIANTS)}")

    config = VARIANTS[variant]
    config["out_png"].parent.mkdir(parents=True, exist_ok=True)
    config["out_csv"].parent.mkdir(parents=True, exist_ok=True)

    base = LocalBase
    base.CENTER_DEG = CENTER_DEG
    model = LocalModel
    sb_df = read_superbubble_table()
    overlay_models, summary = build_overlay_models(variant, model, base, sb_df)
    summary.to_csv(config["out_csv"], index=False, encoding="utf-8-sig")

    fig, axes = plt.subplots(3, 1, figsize=(13.5, 13.1), dpi=360, sharex=True, sharey=True)
    _, mwisp_b_text, mwisp_vmin, mwisp_vmax = draw_mwisp_od_background(axes[0], base)
    _, cfa_b_text, _, _ = draw_lv_background(axes[1], CFA_CACHE, base)
    _, hi_b_text, hi_vmin, hi_vmax = draw_hi_background(axes[2], base)

    draw_overlay_models(axes[0], overlay_models)
    draw_overlay_models(axes[1], overlay_models)
    draw_overlay_models(axes[2], overlay_models)

    style_axis(base, axes[0], r"MWISP $^{12}\mathrm{CO}$ ($|b| \leq 5.25^\circ$)")
    style_axis(base, axes[1], r"CfA $^{12}\mathrm{CO}$ ($|b| \leq 15^\circ$)")
    style_axis(base, axes[2], r"HI4PI H$\,$I ($|b| \leq 15^\circ$)")
    axes[2].set_xlabel("Galactic Longitude l (deg)")

    handles = [
        Line2D(
            [0],
            [0],
            color="0.18",
            lw=0.9,
            linestyle="-",
            label="PV contour from the shell local circular orbit velocity",
        ),
        Line2D(
            [0],
            [0],
            color="0.18",
            lw=0.9,
            linestyle="--",
            dashes=(4, 2),
            label="PV contour calculated from the velocity model",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.018),
        fontsize=LEGEND_FONTSIZE,
        ncol=2,
        frameon=True,
        framealpha=0.9,
    )

    cax_co = fig.add_axes([0.91, 0.42, 0.022, 0.42])
    sm_co = plt.cm.ScalarMappable(norm=LogNorm(vmin=mwisp_vmin, vmax=mwisp_vmax), cmap=plt.get_cmap(BG_CMAP))
    sm_co.set_array([])
    cbar_co = fig.colorbar(sm_co, cax=cax_co)
    cbar_co.set_label(r"$\int T_B\,db$ [K deg]", fontsize=COLORBAR_LABEL_FONTSIZE)
    cbar_co.ax.tick_params(labelsize=COLORBAR_TICK_FONTSIZE)

    cax_hi = fig.add_axes([0.91, 0.13, 0.022, 0.20])
    sm_hi = plt.cm.ScalarMappable(norm=Normalize(vmin=hi_vmin, vmax=hi_vmax, clip=True), cmap=plt.get_cmap(HI_BG_CMAP))
    sm_hi.set_array([])
    cbar_hi = fig.colorbar(sm_hi, cax=cax_hi)
    cbar_hi.set_label(r"$\int T_B\,db$ [K deg]", fontsize=COLORBAR_LABEL_FONTSIZE)
    cbar_hi.ax.tick_params(labelsize=COLORBAR_TICK_FONTSIZE)

    fig.suptitle(config["title"], fontsize=MAIN_TITLE_FONTSIZE, y=0.988)
    fig.subplots_adjust(left=0.075, right=0.89, bottom=0.115, top=0.91, hspace=0.24)
    fig.savefig(config["out_png"], bbox_inches="tight")
    plt.close(fig)
    print(f"Saved figure to: {config['out_png']}")
    print(f"Saved summary to: {config['out_csv']}")
    print(f"Projected {len(overlay_models)} final compromise-ellipse superbubbles with variant={variant}")
    return config["out_png"]


def render_all_variants() -> list[Path]:
    outputs = [render_reference_panel()]
    for variant in VARIANTS:
        outputs.append(render_variant(variant))
    return outputs


def main() -> None:
    """Run Script 21 from validated inputs to the documented outputs."""
    render_all_variants()

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
