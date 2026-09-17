"""Script 22: bubbles gas PV diagram

Purpose
-------
Overlay low-latitude Grade A/B bubbles on MWISP CO PV diagrams.

Method overview
---------------
1. Load the compact MWISP/CfA CO backgrounds and select the configured low-latitude
   Grade A/B bubbles.
2. Use the shared kinematic model to project bubble boundaries into longitude-velocity
   space, splitting wrapped paths at plot boundaries.
3. Render the selected bubble overlays and export the plotted sample information.

Main inputs
-----------
- ../data/COdata/MWISP_12CO_clustered_MCs_PV_projection.fits
- ../data/COdata/CfA_12CO_PV_cache.npz
- ../data/Bubbles.csv

Main outputs
------------
- ../results/figures/22_gradeAB_low_lat_bubbles_PV_overlay.png; closed-bubble PV comparison.

Figure/table role
-----------------
../results/figures/22_gradeAB_low_lat_bubbles_PV_overlay.png; closed-bubble PV comparison.

Runtime and data notes
----------------------
Reference runtime: 5.61 s. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

External-data caveat
--------------------
This script uses packaged intermediate or final inputs unless the inputs section explicitly names an external cache or survey product.

Reading the code
----------------
Start with the path and scientific settings below, then follow render() at
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
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
import matplotlib.patheffects as patheffects
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
os.chdir(HERE)
DATA_DIR = Path("..") / "data"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "22_gradeAB_low_latitude_PV_overlay"
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"
CO_CACHE_DIR = DATA_DIR / "COdata"
OD_FITS = DATA_DIR / "COdata" / "MWISP_12CO_clustered_MCs_PV_projection.fits"
CFA_CACHE = CO_CACHE_DIR / "CfA_12CO_PV_cache.npz"
BUBBLE_CSV = DATA_DIR / "Bubbles.csv"

OUT_PATH = FINAL_FIG_DIR / "22_gradeAB_low_lat_bubbles_PV_overlay.png"
TABLE_PATH = OUT_DIR / "22_gradeAB_low_latitude_bubble_xy_PV_projection_table.csv"

CO_CMAP = "Purples"
CO_VMIN = 1e-2
CO_VMAX = 2e3
EXPANSION_KMS = 10.0
Y_LIMIT = (-50.0, 50.0)
X_LIMIT = (220.0, 10.0)
X_TICKS = np.array([220, 180, 140, 100, 60, 20, 10])
MAX_BACKGROUND_LON_BINS = 7000
PANEL_TITLE_FONTSIZE = 15
MAIN_TITLE_FONTSIZE = 16
LEGEND_FONTSIZE = 14
BUBBLE_LABEL_FONTSIZE = 4.5
COLORBAR_LABEL_FONTSIZE = 14
COLORBAR_TICK_FONTSIZE = 13

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
CENTER_DEG = 180.0


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
    r_hat_y = np.divide(y_gc, r_gal, out=np.zeros_like(x_gc), where=r_gal > 0.0)
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
    reid2019_urc = staticmethod(reid2019_urc)
    lbd_to_vlsr = staticmethod(lbd_to_vlsr)
    plot_longitude_for_center = staticmethod(plot_longitude_for_center)
    absolute_longitude_ticklabels = staticmethod(absolute_longitude_ticklabels)
    split_original_path = staticmethod(split_original_path)


class LocalModel:
    load_base_module = staticmethod(lambda: LocalBase)
    gc_velocity_to_vlsr = staticmethod(gc_velocity_to_vlsr)
    circular_velocity_vector_at_heliocentric_xy = staticmethod(circular_velocity_vector_at_heliocentric_xy)


def load_model_module():
    return LocalModel


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
    vmin = max(float(np.percentile(positive, 2.0)), 1e-3)
    vmax = float(np.percentile(positive, 99.8))
    if not np.isfinite(vmax) or vmax <= vmin:
        vmax = max(vmin * 10.0, float(np.nanmax(positive)))
    return vmin, vmax


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


def draw_od_background(ax, base):
    from astropy.io import fits

    with fits.open(OD_FITS, memmap=True) as hdul:
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

    lon_plot = base.plot_longitude_for_center(lon, center_deg=base.CENTER_DEG)
    order = np.argsort(lon_plot)
    lon_sorted = lon_plot[order]
    lv_sorted = lv_map[:, order]
    lon_sorted, lv_sorted = downsample_longitude(lon_sorted, lv_sorted)

    disp = np.ma.array(np.array(lv_sorted, copy=True), mask=(~np.isfinite(lv_sorted)) | (lv_sorted <= 0.0))
    cmap = plt.get_cmap(CO_CMAP).copy()
    cmap.set_bad((1, 1, 1, 0))
    im = ax.pcolormesh(
        axis_edges(lon_sorted),
        axis_edges(vel),
        disp,
        cmap=cmap,
        norm=LogNorm(vmin=CO_VMIN, vmax=CO_VMAX),
        shading="auto",
        alpha=0.86,
        zorder=0,
    )
    return im


def draw_cfa_background(ax, base):
    with np.load(CFA_CACHE) as cached:
        lv_map = cached["lv_map"].astype(np.float32, copy=False)
        lon = cached["lon"].astype(np.float64, copy=False)
        vel = cached["vel"].astype(np.float64, copy=False)
        b_limit = float(cached["b_limit"][0])

    lon_plot = base.plot_longitude_for_center(lon, center_deg=base.CENTER_DEG)
    order = np.argsort(lon_plot)
    lon_sorted = lon_plot[order]
    lv_sorted = lv_map[:, order]
    disp = np.ma.array(np.array(lv_sorted, copy=True), mask=(~np.isfinite(lv_sorted)) | (lv_sorted <= 0.0))
    vmin, vmax = robust_log_limits(lv_sorted)
    cmap = plt.get_cmap(CO_CMAP).copy()
    cmap.set_bad((1, 1, 1, 0))
    im = ax.pcolormesh(
        axis_edges(lon_sorted),
        axis_edges(vel),
        disp,
        cmap=cmap,
        norm=LogNorm(vmin=vmin, vmax=vmax),
        shading="auto",
        alpha=0.86,
        zorder=0,
    )
    return im, b_limit


def read_grade_ab_bubbles() -> pd.DataFrame:
    df = pd.read_csv(BUBBLE_CSV)
    grade_series = df["Grade"].astype(str).str.upper().str.strip()
    mask = (
        np.isfinite(df["l"])
        & np.isfinite(df["b"])
        & np.isfinite(df["best_distance"])
        & np.isfinite(df["physical_size"])
        & (np.abs(df["b"]) < 5.0)
        & (df["best_distance"] > 0.0)
        & (df["physical_size"] > 0.0)
        & grade_series.isin({"A", "B"})
    )
    if "calc_failed" in df.columns:
        mask &= ~df["calc_failed"].astype(str).str.lower().isin(["true", "1", "yes", "y"])

    out = df[mask].copy().reset_index(drop=True)
    l_rad = np.radians(out["l"].to_numpy(float))
    b_rad = np.radians(out["b"].to_numpy(float))
    d = out["best_distance"].to_numpy(float)
    d_proj = d * np.cos(b_rad)
    out["x_kpc"] = d_proj * np.cos(l_rad)
    out["y_kpc"] = d_proj * np.sin(l_rad)
    out["z_kpc"] = d * np.sin(b_rad)
    out["radius_kpc"] = out["physical_size"].to_numpy(float) / 2000.0
    return out


def bubble_boundary_to_lv(model, base, row, n_samples: int = 181):
    theta = np.linspace(0.0, 2.0 * np.pi, n_samples)
    xc = float(row["x_kpc"])
    yc = float(row["y_kpc"])
    zc = float(row["z_kpc"])
    radius = float(row["radius_kpc"])

    x = xc + radius * np.cos(theta)
    y = yc + radius * np.sin(theta)
    z = np.full_like(x, zc)
    d = np.sqrt(x**2 + y**2 + z**2)
    l_deg = (np.degrees(np.arctan2(y, x)) + 360.0) % 360.0
    b_deg = np.degrees(np.arctan2(z, np.hypot(x, y)))

    v_local = base.lbd_to_vlsr(l_deg, b_deg, d)

    center_vx_gc, center_vy_gc, center_vz_gc = model.circular_velocity_vector_at_heliocentric_xy(base, xc, yc)
    dx = x - xc
    dy = y - yc
    dr = np.hypot(dx, dy)
    ux_h = np.divide(dx, dr, out=np.zeros_like(dx), where=dr > 0.0)
    uy_h = np.divide(dy, dr, out=np.zeros_like(dy), where=dr > 0.0)
    v_expanded = model.gc_velocity_to_vlsr(
        base,
        l_deg,
        b_deg,
        center_vx_gc + EXPANSION_KMS * uy_h,
        center_vy_gc - EXPANSION_KMS * ux_h,
        np.zeros_like(ux_h) + center_vz_gc,
    )

    invalid = (~np.isfinite(d)) | (d <= 0.0) | (~np.isfinite(v_local)) | (~np.isfinite(v_expanded))
    v_local[invalid] = np.nan
    v_expanded[invalid] = np.nan
    return x, y, l_deg, v_local, v_expanded


def in_plot_bounds(l_plot: np.ndarray | float, v_lsr: np.ndarray | float) -> np.ndarray:
    x_min, x_max = sorted(X_LIMIT)
    y_min, y_max = Y_LIMIT
    l_plot = np.asarray(l_plot, dtype=float)
    v_lsr = np.asarray(v_lsr, dtype=float)
    return (
        np.isfinite(l_plot)
        & np.isfinite(v_lsr)
        & (l_plot >= x_min)
        & (l_plot <= x_max)
        & (v_lsr >= y_min)
        & (v_lsr <= y_max)
    )


def visible_path_segments(l_plot: np.ndarray, v_lsr: np.ndarray):
    l_plot = np.asarray(l_plot, dtype=float)
    v_lsr = np.asarray(v_lsr, dtype=float)
    visible = in_plot_bounds(l_plot, v_lsr)
    if not np.any(visible):
        return []
    split_idx = np.where(~visible[:-1] | ~visible[1:])[0] + 1
    out = []
    for ll, vv, ok in zip(np.split(l_plot, split_idx), np.split(v_lsr, split_idx), np.split(visible, split_idx)):
        if ok.size >= 2 and bool(np.all(ok)):
            out.append((ll, vv))
    return out


def plot_lv_bubbles(ax, model, base, bubbles: pd.DataFrame, annotate_ids: bool = True) -> pd.DataFrame:
    summary_rows = []
    for _, row in bubbles.iterrows():
        x, y, l_deg, v_local, v_expanded = bubble_boundary_to_lv(model, base, row)
        for ll, vv in base.split_original_path(x, y, l_deg, v_local):
            for ll_vis, vv_vis in visible_path_segments(ll, vv):
                ax.plot(ll_vis, vv_vis, color="#d62728", lw=0.48, alpha=0.55, linestyle="-", zorder=4, clip_on=True)
        for ll, vv in base.split_original_path(x, y, l_deg, v_expanded):
            for ll_vis, vv_vis in visible_path_segments(ll, vv):
                ax.plot(ll_vis, vv_vis, color="#d62728", lw=0.6, alpha=0.84, linestyle="--", dashes=(3, 2), zorder=5, clip_on=True, path_effects=[patheffects.Stroke(linewidth=0.8, foreground="black", alpha=0.68), patheffects.Normal()])

        center_v = float(
            base.lbd_to_vlsr(
                np.array([float(row["l"])]),
                np.array([float(row["b"])]),
                np.array([float(row["best_distance"])]),
            )[0]
        )
        center_l_plot = float(base.plot_longitude_for_center(np.array([float(row["l"])]), center_deg=base.CENTER_DEG)[0])
        if annotate_ids and bool(in_plot_bounds(center_l_plot, center_v)):
            ax.text(
                center_l_plot,
                center_v,
                str(row["ID"]),
                fontsize=BUBBLE_LABEL_FONTSIZE,
                color="black",
                ha="center",
                va="center",
                alpha=0.82,
                zorder=8,
                path_effects=[patheffects.Stroke(linewidth=0.8, foreground="white", alpha=0.78), patheffects.Normal()],
                clip_on=True,
            )

        summary_rows.append(
            {
                "ID": row["ID"],
                "Grade": row.get("Grade", ""),
                "l_deg": float(row["l"]),
                "b_deg": float(row["b"]),
                "best_distance_kpc": float(row["best_distance"]),
                "physical_size_pc": float(row["physical_size"]),
                "radius_kpc": float(row["radius_kpc"]),
                "x_kpc": float(row["x_kpc"]),
                "y_kpc": float(row["y_kpc"]),
                "z_kpc": float(row["z_kpc"]),
                "v_center_model_kms": center_v,
                "v_local_min_kms": float(np.nanmin(v_local)),
                "v_local_max_kms": float(np.nanmax(v_local)),
                "v_expanded_min_kms": float(np.nanmin(v_expanded)),
                "v_expanded_max_kms": float(np.nanmax(v_expanded)),
            }
        )
    return pd.DataFrame(summary_rows)


def style_lv_axis(base, ax, title: str):
    ax.axhline(0.0, color="0.45", lw=0.7, alpha=0.5, zorder=1)
    ax.set_xlim(*X_LIMIT)
    ax.set_ylim(-50.0, 50.0)
    ax.set_xticks(X_TICKS)
    ax.set_xticklabels([str(int(tick)) for tick in X_TICKS])
    ax.grid(color="0.86", linestyle="--", linewidth=0.45, alpha=0.9)
    ax.set_title(title, fontsize=PANEL_TITLE_FONTSIZE)
    ax.set_ylabel(r"$V_{\rm LSR}$ (km s$^{-1}$)")


def render() -> tuple[Path, Path]:
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    TABLE_PATH.parent.mkdir(parents=True, exist_ok=True)

    model = load_model_module()
    base = model.load_base_module()
    bubbles = read_grade_ab_bubbles()

    fig, ax = plt.subplots(1, 1, figsize=(13.5, 6.0), dpi=450)
    draw_od_background(ax, base)

    summary = plot_lv_bubbles(ax, model, base, bubbles, annotate_ids=True)

    title_suffix = "Grade A+B"
    style_lv_axis(
        base,
        ax,
        rf"MWISP $^{{12}}\mathrm{{CO}}$ ($|b| \leq 5.25^\circ$), {title_suffix} bubbles (N={len(bubbles)}, $|b| < 5^\circ$)",
    )
    ax.set_xlabel("Galactic Longitude l (deg)")

    handles = [
        Line2D(
            [0],
            [0],
            color="#d62728",
            lw=1.0,
            linestyle="-",
            label="PV contour from the shell local circular orbit velocity",
        ),
        Line2D(
            [0],
            [0],
            color="#d62728",
            lw=2.4,
            linestyle="--",
            dashes=(3, 2),
            label="PV contour calculated from the velocity model",
            path_effects=[patheffects.Stroke(linewidth=0.8, foreground="black", alpha=0.68), patheffects.Normal()],
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

    fig.subplots_adjust(left=0.075, right=0.91, bottom=0.26, top=0.90)
    cax_bg = fig.add_axes([0.935, 0.32, 0.022, 0.50])
    sm_bg = plt.cm.ScalarMappable(norm=LogNorm(vmin=CO_VMIN, vmax=CO_VMAX), cmap=plt.get_cmap(CO_CMAP))
    sm_bg.set_array([])
    cbar_bg = fig.colorbar(sm_bg, cax=cax_bg)
    cbar_bg.set_label(r"$\int T_B\,db$ [K deg]", fontsize=COLORBAR_LABEL_FONTSIZE)
    cbar_bg.ax.tick_params(labelsize=COLORBAR_TICK_FONTSIZE)


    fig.savefig(OUT_PATH, bbox_inches="tight")
    plt.close(fig)
    try:
        summary.to_csv(TABLE_PATH, index=False, encoding="utf-8-sig")
    except PermissionError as exc:
        print(f"Warning: could not overwrite projection table: {exc}")

    print(f"Selected {len(bubbles)} Grade A+B bubbles with |b| < 5 deg")
    print(f"Saved PV figure to: {OUT_PATH}")
    print(f"Saved projection table to: {TABLE_PATH}")
    return OUT_PATH, TABLE_PATH

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
    _run_with_timing(render, __file__)
