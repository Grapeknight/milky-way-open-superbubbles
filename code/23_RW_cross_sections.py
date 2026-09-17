#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Render the Radcliffe-Wave dust slice and young-cluster kinematics as Fig. 3.

The slice is a vertical plane whose XY trace is y=tan(angle)*x+intercept.
At the defaults (60 degrees, intercept 0.6 kpc), a band of half-width
0.1 kpc supplies the dust and cluster samples. T increases along
(cos(angle), sin(angle)); its zero is the perpendicular foot from the XY
origin to the trace. The normal n=(-sin(angle), cos(angle)) sets signed
slice distances and transverse velocities. Z retains the input vertical
coordinate. Geometry and dust selection use kpc, plotted axes use pc,
and the cluster table supplies positions in pc and LSR velocities in km/s.

Dust comes from the smoothed 3D parquet cube, falling back to the raw cube.
For the default slice, the CSV cache under
../results/intermediate_output/23_radcliffe_wave_fig3_slices/ is reused or
created when missing; changed slice parameters are computed from the cube.
In each T-Z bin, the mean dust density is multiplied by the full band
thickness (2*half_width) to approximate Delta E(B-V) in mag. This is a
slab-integrated map of a finite band, whereas the shell outlines are
intersections with its zero-width central plane.

Shell geometry is read from ../results/superbubble_final_fit_parameters.csv.
Rotated ellipses yield analytic T-Z ellipsoid intersections; cylinders use
a fixed 0.05 kpc half-height. SB21 uses its tabulated cap-plane geometry.
All other ellipsoids, including SB31, retain their original fitted axes and
boundaries with no display-specific radius adjustment. Marks 1 and 2 show
the lower and upper arcs, respectively; other marks show full ellipses.
Solid/dashed outlines identify centers at non-positive/positive signed
normal distance, respectively (the figure's front/behind convention).

The cluster input is Script 13's Hunt_Konietzka_master_table_compact_columns.csv.
Complete rows younger than 30 Myr inside the same band and zoom field are
projected onto T and n. Marker area increases with age; arrows show (v_T,v_Z),
and marker/arrow colors encode v_n on a common -20 to +20 km/s scale.
Konietzka and Hunt sources use stars and circles. The fixed Orion-Eridanus
circle is a visual reference, not an intersection derived from the fit table.

../results/figures/23_fig3.png is written directly as one combined figure:
(a) a dust overview and (b) a grayscale dust/cluster zoom. A dashed overview
box and connectors identify the zoom boundaries. The Radcliffe Wave label
is drawn above these overlays. SB17 is omitted from the lower panel only.
All paths are relative to code/; the upstream fits and cluster data are read
without modification."""

from __future__ import annotations
import time

import argparse
import os
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.ticker import FormatStrFormatter
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, ConnectionPatch, Rectangle
import numpy as np
import pandas as pd
import pyarrow.parquet as pq


# ----------------------------------------------------------------------------
# Path configuration; paths are relative to code/ and point directly into ../data or ../results/.
# ----------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
os.chdir(HERE)
FINAL_DIR = Path("..") / "results"
FINAL_FIG_DIR = FINAL_DIR / "figures"
OUT_DIR = Path("..") / "results" / "intermediate_output" / "23_radcliffe_wave_fig3_slices"
DUST_DATA_DIR = Path("..") / "results" / "intermediate_output" / "3d_dust_map_products"
DATA_DIR = Path("..") / "data"

DUST_PARQUET = DUST_DATA_DIR / "smoothed_3d_dust_cube.parquet"
FALLBACK_DUST_PARQUET = DUST_DATA_DIR / "raw_3d_dust_cube.parquet"
DUST_S_CSV = OUT_DIR / "radcliffe_wave_dust_slice_angle60_b0p60_hw0p10.csv"
KONIETZKA_CSV = DATA_DIR / "star_cluster_data" / "Konietzka2023.csv"
HUNT_KONIETZKA_MINIMAL_CSV = (
    Path("..") / "results" / "intermediate_output"
    / "13_traceback_cluster_input_data"
    / "Hunt_Konietzka_master_table_compact_columns.csv"
)
FINAL_TABLE = FINAL_DIR / "superbubble_final_fit_parameters.csv"

OUT_FIG = FINAL_FIG_DIR / "23_fig3.png"


# ----------------------------------------------------------------------------
# Slice definition in the input XY frame; plotted T and Z are converted to pc.
# ----------------------------------------------------------------------------
ANGLE_DEG = 60.0
INTERCEPT_KPC = 0.6
HALF_WIDTH_KPC = 0.1

T_BINS = 301
Z_BINS = 101
TZ_TOP_XLIM_PC = (-3500, 3300)
TZ_TOP_YLIM_PC = (-500, 500)
TZ_BOTTOM_XLIM_PC = (-1500, 1000)
TZ_BOTTOM_YLIM_PC = (-400, 400)
BOTTOM_CLUSTER_AGE_CUT_MYR = 30.0

DUST_FILTER_XY_KPC = (-3.6, 3.6)
DUST_FILTER_Z_KPC = (-0.6, 0.6)
DUST_DENSITY_VMIN = 0.0
DUST_DENSITY_VMAX = 0.8
DELTA_EBV_VMIN = 0.0
DELTA_EBV_VMAX = 0.2

CYLINDER_HALF_HEIGHT_KPC = 0.05
EPS = 1e-12
MARK_COLORS = {1: "blue", 2: "red", 3: "green"}
MARK_LABELS = {
    1: "N-Open",
    2: "S-Open",
    3: "Disk-penetrating",
}
SB21_CAP_SECTION_IDS = {21}

AXIS_LABEL_FONTSIZE = 13
TICK_LABEL_FONTSIZE = 10
COLORBAR_LABEL_FONTSIZE = 11
COLORBAR_TICK_FONTSIZE = 10
RW_LABEL_FONTSIZE = 11
SB_LABEL_FONTSIZE = 7.5
PANEL_LABEL_FONTSIZE = 11
ORION_LABEL_FONTSIZE = 7.5
ORION_COLOR = "darkviolet"
ORION_LINESTYLE = "solid"
ORION_DRAW_DISTANCE_KPC = -np.inf
ORION_CENTER_PC = (-237, -170)
ORION_RADIUS_PC = 80
LEGEND_FONTSIZE = 7.5
BOTTOM_LEGEND_FONTSIZE = 8.5
FRONT_LINESTYLE = "solid"
BACK_LINESTYLE = (0, (4.0, 3.0))
SB_LABEL_Z_OFFSETS_PC = {
    15: -35.0,
    16: 35.0,
    17: -35.0,
    30: -35.0,
}


def parse_args() -> argparse.Namespace:
    """Read slice angle (degrees), Y intercept (kpc), band half-width (kpc), and PNG dpi."""
    parser = argparse.ArgumentParser(description='Run Script 23: RW cross sections.')
    parser.add_argument("--angle", type=float, default=ANGLE_DEG, help='Command-line option for the documented workflow.')
    parser.add_argument("--intercept", type=float, default=INTERCEPT_KPC, help='Command-line option for the documented workflow.')
    parser.add_argument("--hw", type=float, default=HALF_WIDTH_KPC, help='Command-line option for the documented workflow.')
    parser.add_argument("--dpi", type=int, default=300, help='Command-line option for the documented workflow.')
    return parser.parse_args()


def mark_color(mark) -> str:
    """Use the common black shell color; front/behind is encoded by line style."""
    return "black"


def line_geometry(angle_deg: float, intercept_kpc: float):
    """Return angle, tangent/normal unit vectors, and the T=0 point in kpc.

    The closest point to the origin is intercept*cos(theta)*unit_n. Hence
    T=dot(XY-point_t0, unit_t), while signed normal distance is
    -dot(XY, (sin(theta), -cos(theta)))-intercept*cos(theta)."""
    theta = np.radians(angle_deg)
    unit_t = np.array([np.cos(theta), np.sin(theta)])
    unit_n = np.array([-np.sin(theta), np.cos(theta)])
    point_t0 = np.array([
        -intercept_kpc * np.sin(theta) * np.cos(theta),
        intercept_kpc * np.cos(theta) ** 2,
    ])
    return theta, unit_t, unit_n, point_t0


def read_xyz_dust_subset(parquet_path: Path) -> pd.DataFrame:
    """Load only the local XYZ box and normalize supported parquet column names.

    Coordinates are in kpc. The dust/dust_raw column is the differential
    reddening density used below in mag/kpc; parquet filtering limits the
    read to the configured XY and Z bounds before slicing."""
    names = set(pq.ParquetFile(str(parquet_path)).schema.names)
    if {"x", "y", "z", "dust"}.issubset(names):
        cols = ["x", "y", "z", "dust"]
        rename = {}
    elif {"X", "Y", "Z", "dust"}.issubset(names):
        cols = ["X", "Y", "Z", "dust"]
        rename = {"X": "x", "Y": "y", "Z": "z"}
    elif {"x", "y", "z", "dust_raw"}.issubset(names):
        cols = ["x", "y", "z", "dust_raw"]
        rename = {"dust_raw": "dust"}
    elif {"X", "Y", "Z", "dust_raw"}.issubset(names):
        cols = ["X", "Y", "Z", "dust_raw"]
        rename = {"X": "x", "Y": "y", "Z": "z", "dust_raw": "dust"}
    else:
        raise ValueError(f"Could not identify XYZ parquet columns: {sorted(names)}")
    x_col, y_col, z_col = cols[0], cols[1], cols[2]
    return pd.read_parquet(
        parquet_path,
        columns=cols,
        filters=[
            (x_col, ">=", DUST_FILTER_XY_KPC[0]), (x_col, "<=", DUST_FILTER_XY_KPC[1]),
            (y_col, ">=", DUST_FILTER_XY_KPC[0]), (y_col, "<=", DUST_FILTER_XY_KPC[1]),
            (z_col, ">=", DUST_FILTER_Z_KPC[0]), (z_col, "<=", DUST_FILTER_Z_KPC[1]),
        ],
    ).rename(columns=rename)


def load_dust_subset() -> pd.DataFrame:
    """Prefer the smoothed cube; use the raw cube when the smoothed file is absent."""
    if DUST_PARQUET.exists():
        return read_xyz_dust_subset(DUST_PARQUET)
    return read_xyz_dust_subset(FALLBACK_DUST_PARQUET)


def build_band_slice_df(
    df: pd.DataFrame,
    theta: float,
    point_t0: np.ndarray,
    intercept_kpc: float,
    half_width: float,
) -> pd.DataFrame:
    """Select |normal distance| <= half_width and project each cell along T.

    Keep the original XYZ coordinates alongside dust and the projected s
    column (the plotted T coordinate). All distances here remain in kpc."""
    dist = (
        -np.sin(theta) * df["x"].to_numpy()
        + np.cos(theta) * df["y"].to_numpy()
        - intercept_kpc * np.cos(theta)
    )
    mask = np.abs(dist) <= half_width
    cut = df.loc[mask, ["x", "y", "z", "dust"]]
    if cut.empty:
        raise ValueError("No dust points are available in the current slice; check the angle, intercept, or half-width.")

    result = cut.rename(columns={"x": "X", "y": "Y", "z": "Z"}).copy()
    result["s"] = (
        (result["X"].to_numpy() - point_t0[0]) * np.cos(theta)
        + (result["Y"].to_numpy() - point_t0[1]) * np.sin(theta)
    )
    return result[["X", "Y", "Z", "dust", "s"]]


def load_or_create_rw_dust_slice(df: pd.DataFrame | None = None) -> pd.DataFrame:
    """Reuse the default 60-degree slice CSV or build it from the available cube.

    The filename describes only the fixed default geometry. main() bypasses
    this cache whenever angle, intercept, or half-width differs from defaults."""
    if DUST_S_CSV.exists():
        return pd.read_csv(DUST_S_CSV).dropna(subset=["dust", "s", "Z"])

    if df is None:
        df = read_xyz_dust_subset(DUST_PARQUET if DUST_PARQUET.exists() else FALLBACK_DUST_PARQUET)
    theta, _, _, point_t0 = line_geometry(ANGLE_DEG, INTERCEPT_KPC)
    band = build_band_slice_df(df, theta, point_t0, INTERCEPT_KPC, HALF_WIDTH_KPC)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    band.to_csv(DUST_S_CSV, index=False, float_format="%.8g")
    print(f"Saved RW dust slice cache: {DUST_S_CSV}")
    return band


def build_dust_tz_map_from_slice(slice_df: pd.DataFrame):
    """Bin finite dust densities into T-Z cells and return their arithmetic means.

    T_BINS and Z_BINS count edges, so the default image has 300 by 100 cells.
    Empty bins remain NaN. Transpose the histogram to Z rows and T columns
    for pcolormesh; conversion from density to reddening occurs in main()."""
    s = slice_df["s"].to_numpy()
    z = slice_df["Z"].to_numpy()
    dust = slice_df["dust"].to_numpy()
    valid = np.isfinite(dust)

    t_edges = np.linspace(np.nanmin(s[valid]), np.nanmax(s[valid]), T_BINS)
    z_edges = np.linspace(np.nanmin(z[valid]), np.nanmax(z[valid]), Z_BINS)
    hist, _, _ = np.histogram2d(s[valid], z[valid], bins=[t_edges, z_edges], weights=dust[valid])
    counts, _, _ = np.histogram2d(s[valid], z[valid], bins=[t_edges, z_edges])
    mean_dust = np.divide(
        hist,
        counts,
        out=np.full_like(hist, np.nan, dtype=float),
        where=counts != 0,
    ).T
    return t_edges, z_edges, mean_dust


def build_dust_tz_map(
    df: pd.DataFrame,
    theta: float,
    point_t0: np.ndarray,
    intercept_kpc: float,
    half_width: float,
    t_range_kpc: tuple[float, float] | None = None,
    z_range_kpc: tuple[float, float] | None = None,
):
    """Build a band and bin its full occupied extent.

    The optional t_range_kpc/z_range_kpc arguments are retained by this helper
    but are not applied; bin edges follow the selected data extent."""
    band = build_band_slice_df(df, theta, point_t0, intercept_kpc, half_width)
    return build_dust_tz_map_from_slice(band)


def load_superbubbles() -> pd.DataFrame:
    """Read fitted shell and XY-plane geometry; coerce numeric fields without refitting."""
    needed = {
        "id", "mark", "shape",
        "center_x_kpc", "center_y_kpc", "center_z_kpc",
        "a_radius_kpc", "b_radius_kpc", "c_radius_kpc", "angle_deg",
        "xy_plane_center_x_kpc", "xy_plane_center_y_kpc",
        "xy_plane_z_kpc", "xy_plane_a_kpc", "xy_plane_b_kpc", "xy_plane_angle_deg",
    }
    df = pd.read_csv(FINAL_TABLE, usecols=lambda c: c in needed)
    for col in df.columns:
        if col != "shape":
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def xy_terms(center_xy: np.ndarray, axis_a: float, axis_b: float, angle_deg: float, unit_t, point_t0):
    """Substitute XY=point_t0+T*unit_t into the rotated ellipse equation.

    Return A, B, C for A*T**2+2*B*T+C, the normalized squared XY distance
    from the shell center. Axis lengths and positions must use the same
    kpc units. A nonfinite angle falls back to an unrotated ellipse."""
    if not np.isfinite(axis_a) or not np.isfinite(axis_b) or axis_a <= 0 or axis_b <= 0:
        return None
    phi = np.radians(float(angle_deg) if np.isfinite(angle_deg) else 0.0)
    major = np.array([np.cos(phi), np.sin(phi)])
    minor = np.array([-np.sin(phi), np.cos(phi)])
    delta = point_t0 - center_xy

    alpha_maj = float(np.dot(unit_t, major))
    alpha_min = float(np.dot(unit_t, minor))
    beta_maj = float(np.dot(delta, major))
    beta_min = float(np.dot(delta, minor))

    qa = 1.0 / (axis_a * axis_a)
    qb = 1.0 / (axis_b * axis_b)
    A = alpha_maj ** 2 * qa + alpha_min ** 2 * qb
    B = alpha_maj * beta_maj * qa + alpha_min * beta_min * qb
    C = beta_maj ** 2 * qa + beta_min ** 2 * qb
    if A <= EPS:
        return None
    return A, B, C


def xy_min_from_terms(terms):
    """Return C-B**2/A, the ellipse quadratic minimum at T=-B/A."""
    A, B, C = terms
    return C - B * B / A


def choose_slice_position(center_xy, axis_a, axis_b, angle, unit_t, unit_n, point_t0, half_width):
    """Prefer the central plane; optionally find a plane within a finite slab.

    If only the slab intersects, minimize the XY quadratic over normal
    offset using endpoints and a bounded ternary search. main() passes
    half_width=0, so this optional branch is not used by the final figure."""
    def terms_at(offset):
        shifted_point_t0 = point_t0 + offset * unit_n
        terms = xy_terms(center_xy, axis_a, axis_b, angle, unit_t, shifted_point_t0)
        if terms is None:
            return None, np.inf
        return terms, xy_min_from_terms(terms)

    central_terms, central_xy_min = terms_at(0.0)
    if central_terms is not None and central_xy_min <= 1.0 + 1e-9:
        return central_terms, 0.0, central_xy_min
    if half_width <= EPS:
        return central_terms, 0.0, central_xy_min

    candidates = []
    for offset in (-half_width, half_width):
        terms, xy_min = terms_at(offset)
        candidates.append((xy_min, offset, terms))

    lo, hi = -half_width, half_width
    for _ in range(64):
        left = lo + (hi - lo) / 3.0
        right = hi - (hi - lo) / 3.0
        _, left_xy_min = terms_at(left)
        _, right_xy_min = terms_at(right)
        if left_xy_min < right_xy_min:
            hi = right
        else:
            lo = left
    offset = 0.5 * (lo + hi)
    terms, xy_min = terms_at(offset)
    candidates.append((xy_min, offset, terms))

    best_xy_min, best_offset, best_terms = min(candidates, key=lambda item: item[0])
    return best_terms, best_offset, best_xy_min


def compute_cut(
    row: pd.Series,
    unit_t: np.ndarray,
    point_t0: np.ndarray,
    unit_n: np.ndarray | None = None,
    half_width: float = 0.0,
    extra_cap_section_ids: set[int] | None = None,
):
    """Return analytic section bounds in kpc, or None when no finite cut exists.

    After completing the XY square, remaining=1-(C-B**2/A) determines the
    available ellipsoid section: its T radius is sqrt(remaining/A) and its
    Z radius is c*sqrt(remaining). Cylinder height is independent of this
    factor. SB21's special cap uses the XY-plane ellipse and base-to-apex
    height; other ellipsoids use their original full fitted axes, including
    SB31. Arc selection by mark is deferred to the drawing helper."""
    shape = str(row["shape"]).strip().lower()
    sb_id = int(row["id"])
    mark = int(row["mark"]) if np.isfinite(row["mark"]) else -1
    cap_section_ids = SB21_CAP_SECTION_IDS | (extra_cap_section_ids or set())

    if shape == "cylinder":
        center_xy = np.array([row["center_x_kpc"], row["center_y_kpc"]], dtype=float)
        axis_a = float(row["a_radius_kpc"])
        axis_b = float(row["b_radius_kpc"])
        angle = float(row["angle_deg"]) if np.isfinite(row["angle_deg"]) else 0.0
        z_center = float(row["center_z_kpc"])
        z_radius = CYLINDER_HALF_HEIGHT_KPC
        cut_kind = "rectangle"
        section_model = "cylinder"
    elif sb_id in cap_section_ids:
        center_xy = np.array([
            row.get("xy_plane_center_x_kpc", row["center_x_kpc"]),
            row.get("xy_plane_center_y_kpc", row["center_y_kpc"]),
        ], dtype=float)
        axis_a = float(row["xy_plane_a_kpc"])
        axis_b = float(row["xy_plane_b_kpc"])
        angle = float(row["xy_plane_angle_deg"]) if np.isfinite(row["xy_plane_angle_deg"]) else float(row["angle_deg"])
        xy_z = float(row["xy_plane_z_kpc"])
        c_full = float(row["c_radius_kpc"])
        z_apex = float(row["center_z_kpc"]) - c_full if mark == 1 else float(row["center_z_kpc"]) + c_full
        z_center = xy_z
        z_radius = abs(xy_z - z_apex)
        cut_kind = "ellipse_arc"
        section_model = "ellipsoid_cap"
    else:
        center_xy = np.array([row["center_x_kpc"], row["center_y_kpc"]], dtype=float)
        axis_a = float(row["a_radius_kpc"])
        axis_b = float(row["b_radius_kpc"])
        angle = float(row["angle_deg"]) if np.isfinite(row["angle_deg"]) else 0.0
        z_center = float(row["center_z_kpc"])
        z_radius = float(row["c_radius_kpc"])
        cut_kind = "ellipse_arc"
        section_model = "full_ellipsoid"

    if unit_n is None:
        unit_n = np.array([-unit_t[1], unit_t[0]], dtype=float)
    terms, slice_offset, xy_min = choose_slice_position(
        center_xy, axis_a, axis_b, angle, unit_t, unit_n, point_t0, half_width
    )
    if terms is None:
        return None
    A, B, C = terms
    t_center = -B / A
    if xy_min > 1.0 + 1e-9:
        return None
    # Clamp roundoff at a tangent; reject zero-size cuts below.
    remaining = max(0.0, 1.0 - xy_min)
    t_radius = float(np.sqrt(remaining / A))
    if cut_kind == "ellipse_arc":
        z_radius = z_radius * np.sqrt(remaining)
    if t_radius <= EPS or z_radius <= EPS:
        return None

    result = {
        "id": sb_id,
        "mark": mark,
        "cut_kind": cut_kind,
        "section_model": section_model,
        "t_center_kpc": t_center,
        "z_center_kpc": z_center,
        "t_radius_kpc": t_radius,
        "z_radius_kpc": z_radius,
        "t_min_kpc": t_center - t_radius,
        "t_max_kpc": t_center + t_radius,
        "z_min_kpc": z_center - z_radius,
        "z_max_kpc": z_center + z_radius,
        "slice_offset_kpc": slice_offset,
        "center_signed_distance_kpc": float(np.dot(center_xy - point_t0, unit_n)),
    }
    if section_model == "full_ellipsoid":
        result.update({
            "full_t_center_kpc": t_center,
            "full_z_center_kpc": z_center,
            "full_t_radius_kpc": t_radius,
            "full_z_radius_kpc": z_radius,
        })
    return result


def compute_all_cuts(
    superbubbles: pd.DataFrame,
    unit_t: np.ndarray,
    point_t0: np.ndarray,
    unit_n: np.ndarray | None = None,
    half_width: float = 0.0,
    extra_cap_section_ids: set[int] | None = None,
) -> pd.DataFrame:
    """Collect intersecting shells, retaining IDs, section models, and signed distances."""
    rows = []
    for _, row in superbubbles.iterrows():
        cut = compute_cut(
            row, unit_t, point_t0, unit_n=unit_n, half_width=half_width,
            extra_cap_section_ids=extra_cap_section_ids,
        )
        if cut is not None:
            rows.append(cut)
    return pd.DataFrame(rows)


def cut_linestyle(cut):
    """Map center distance <= 0 to solid (front), positive to dashed (behind)."""
    signed_distance = cut.get("center_signed_distance_kpc", 0.0)
    return FRONT_LINESTYLE if signed_distance <= 0 else BACK_LINESTYLE


def draw_ellipsoid_cut_pc(ax, cut, color, outline=False):
    """Convert the section to pc and draw its retained cap or full ellipse.

    Mark 1 uses angles pi..2*pi (lower Z), mark 2 uses 0..pi (upper Z).
    The optional solid white underlay improves visibility against dust."""
    mark = int(cut["mark"])
    if mark == 1:
        angles = np.linspace(np.pi, 2 * np.pi, 180)
    elif mark == 2:
        angles = np.linspace(0, np.pi, 180)
    else:
        angles = np.linspace(0, 2 * np.pi, 240)
    t = 1000.0 * (cut["t_center_kpc"] + cut["t_radius_kpc"] * np.cos(angles))
    z = 1000.0 * (cut["z_center_kpc"] + cut["z_radius_kpc"] * np.sin(angles))
    linestyle = cut_linestyle(cut)
    if outline:
        ax.plot(t, z, color="white", linewidth=3.4, alpha=0.95, zorder=6, linestyle=FRONT_LINESTYLE)
    ax.plot(t, z, color=color, linewidth=1.55, alpha=0.98, zorder=7, linestyle=linestyle)


def draw_cylinder_cut_pc(ax, cut, color, outline=False):
    """Render the cylinder section as a rectangle in pc with the same depth convention."""
    x0 = 1000.0 * cut["t_min_kpc"]
    y0 = 1000.0 * cut["z_min_kpc"]
    width = 2000.0 * cut["t_radius_kpc"]
    height = 2000.0 * cut["z_radius_kpc"]
    linestyle = cut_linestyle(cut)
    if outline:
        ax.add_patch(Rectangle((x0, y0), width, height, edgecolor="white", facecolor="none",
                               linewidth=3.4, alpha=0.95, zorder=6, linestyle=FRONT_LINESTYLE))
    ax.add_patch(Rectangle((x0, y0), width, height, edgecolor=color, facecolor="none",
                           linewidth=1.55, alpha=0.98, zorder=7, linestyle=linestyle))


def draw_superbubble_cuts_pc(ax, cuts: pd.DataFrame, xlim, ylim, outline=False):
    """Draw farther/positive-distance centers first, then label visible sections.

    Plot bounds cull off-panel shells and labels. Fixed label offsets affect
    only annotation positions; none of them changes a shell boundary."""
    if "center_signed_distance_kpc" in cuts.columns:
        draw_cuts = cuts.sort_values("center_signed_distance_kpc", ascending=False)
    else:
        draw_cuts = cuts
    for _, cut in draw_cuts.iterrows():
        if 1000.0 * cut["t_max_kpc"] < xlim[0] or 1000.0 * cut["t_min_kpc"] > xlim[1]:
            continue
        if 1000.0 * cut["z_max_kpc"] < ylim[0] or 1000.0 * cut["z_min_kpc"] > ylim[1]:
            continue
        color = mark_color(cut["mark"])
        if cut["cut_kind"] == "rectangle":
            draw_cylinder_cut_pc(ax, cut, color, outline=outline)
        else:
            draw_ellipsoid_cut_pc(ax, cut, color, outline=outline)

        label_z = cut["z_center_kpc"]
        mark = int(cut["mark"])
        if mark == 1:
            label_z -= 0.35 * cut["z_radius_kpc"]
        elif mark == 2:
            label_z += 0.35 * cut["z_radius_kpc"]
        label_z += SB_LABEL_Z_OFFSETS_PC.get(int(cut["id"]), 0.0) / 1000.0
        label_t_pc = 1000.0 * cut["t_center_kpc"]
        label_z_pc = 1000.0 * label_z
        if xlim[0] <= label_t_pc <= xlim[1] and ylim[0] <= label_z_pc <= ylim[1]:
            ax.text(
                label_t_pc, label_z_pc, f"SB{int(cut['id'])}",
                fontsize=SB_LABEL_FONTSIZE, ha="center", va="center", color="black", weight="bold",
                zorder=8, clip_on=True,
                path_effects=[pe.withStroke(linewidth=2.8, foreground="white")],
            )


def draw_orion_eridanus_marker(ax):
    """Draw the fixed reference circle at (T,Z)=(-237,-170) pc with radius 80 pc."""
    center = ORION_CENTER_PC
    radius = ORION_RADIUS_PC
    ax.add_patch(Circle(center, radius=radius, edgecolor="white", facecolor="none",
                        linewidth=3.4, alpha=0.95, zorder=6, linestyle=ORION_LINESTYLE))
    ax.add_patch(Circle(center, radius=radius, edgecolor=ORION_COLOR, facecolor="none",
                        linewidth=1.7, alpha=0.98, zorder=7, linestyle=ORION_LINESTYLE))


def draw_top_panel(ax, t_edges, z_edges, delta_ebv_map, cuts, title: str, norm, cmap: str):
    """Draw the dust overview, fitted shell sections, and Radcliffe Wave span.

    The span equals the lower panel's T limits. Its text uses zorder 40 so
    the later zoom box (19) and connectors (18) cannot cover the label."""
    s_edges_pc = 1000.0 * t_edges
    z_edges_pc = 1000.0 * z_edges

    mesh = ax.pcolormesh(s_edges_pc, z_edges_pc, delta_ebv_map, cmap=cmap,
                         norm=norm, shading="auto", rasterized=True)

    line_z = 320
    rw_t_min, rw_t_max = TZ_BOTTOM_XLIM_PC
    ax.annotate("", xy=(rw_t_min, line_z), xytext=(rw_t_max, line_z),
                arrowprops={"arrowstyle": "<->", "color": "black", "lw": 1.6}, zorder=20)
    line_length = 100
    ax.plot([rw_t_min, rw_t_min], [line_z - line_length / 2, line_z + line_length / 2], "k-", lw=1.6)
    ax.plot([rw_t_max, rw_t_max], [line_z - line_length / 2, line_z + line_length / 2], "k-", lw=1.6)
    ax.text((rw_t_min + rw_t_max) / 2, line_z + line_length / 2 + 10, "Radcliffe Wave",
            ha="center", va="bottom", fontsize=RW_LABEL_FONTSIZE, fontweight="bold", color="black",
            zorder=40, path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])

    ax.set_aspect("equal", adjustable="box")
    draw_superbubble_cuts_pc(ax, cuts, TZ_TOP_XLIM_PC, TZ_TOP_YLIM_PC, outline=True)
    draw_orion_eridanus_marker(ax)

    if title:
        ax.text(0.01, 0.94, title, transform=ax.transAxes, ha="left", va="top",
                fontsize=PANEL_LABEL_FONTSIZE, fontweight="bold", color="black",
                path_effects=[pe.withStroke(linewidth=3.0, foreground="white")])
    ax.set_xlim(*TZ_TOP_XLIM_PC)
    ax.set_ylim(*TZ_TOP_YLIM_PC)
    ax.set_ylabel("Z (pc)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True,
                   labelsize=TICK_LABEL_FONTSIZE)
    ax.grid(True, linestyle="--", alpha=0.5)
    return mesh


def apply_legend_line_effects(legend):
    for line in legend.get_lines():
        line.set_path_effects([pe.withStroke(linewidth=3.2, foreground="white")])


def add_superbubble_legend(ax, *, split_right: bool = False, fontsize: float = LEGEND_FONTSIZE):
    """Explain the front/behind styles and Orion reference, optionally in two legends."""
    black_handles = [
        Line2D([0], [0], color="black", lw=1.8, linestyle=FRONT_LINESTYLE,
               label="Open superbubbles (centers in front of slice)"),
        Line2D([0], [0], color="black", lw=1.8, linestyle=BACK_LINESTYLE,
               label="Open superbubbles (centers behind slice)"),
    ]
    orion_handles = [
        Line2D([0], [0], color=ORION_COLOR, lw=1.8, linestyle=ORION_LINESTYLE,
               label="Orion–Eridanus superbubble"),
    ]
    if split_right:
        black_legend = ax.legend(handles=black_handles, loc="upper right", fontsize=fontsize,
                                 framealpha=0.82, borderpad=0.35, handlelength=2.4,
                                 labelspacing=0.32)
        apply_legend_line_effects(black_legend)
        ax.add_artist(black_legend)
        orion_legend = ax.legend(handles=orion_handles, loc="lower right", fontsize=fontsize,
                                 framealpha=0.82, borderpad=0.35, handlelength=2.4,
                                 labelspacing=0.32)
        apply_legend_line_effects(orion_legend)
        return

    legend = ax.legend(handles=black_handles + orion_handles, loc="upper left",
                       fontsize=fontsize, framealpha=0.82, borderpad=0.35,
                       handlelength=2.4, labelspacing=0.32)
    apply_legend_line_effects(legend)


def draw_figure(t_edges, z_edges, delta_ebv_map, cuts, clusters, dpi: int, ebv_vmin: float, ebv_vmax: float):
    """Assemble overview and cluster zoom into the single manuscript figure.

    The box uses the exact lower-panel T/Z limits in overview data units.
    Connectors join its lower corners to the zoom's top corners using axes
    fractions at the destination, preserving the link under layout changes."""
    fig = plt.figure(figsize=(12, 6.6))
    dust_cbar_ax = fig.add_axes([0.18, 0.91, 0.62, 0.018])
    ax_top = fig.add_axes([0.08, 0.64, 0.86, 0.26])
    ax_bottom = fig.add_axes([0.18, 0.15, 0.66, 0.384])
    velocity_cbar_ax = fig.add_axes([0.24, 0.06, 0.54, 0.018])

    norm = Normalize(vmin=ebv_vmin, vmax=ebv_vmax, clip=True)
    mesh = draw_top_panel(ax_top, t_edges, z_edges, delta_ebv_map, cuts, "(a)",
                          norm, "Spectral_r")
    add_superbubble_legend(ax_top, split_right=True, fontsize=BOTTOM_LEGEND_FONTSIZE)
    ax_top.set_xlabel("T (pc)", fontsize=AXIS_LABEL_FONTSIZE)

    cbar = fig.colorbar(mesh, cax=dust_cbar_ax, orientation="horizontal", extend="max")
    cbar.ax.text(1.08, 0.5, "ΔE(B-V) [mag]", transform=cbar.ax.transAxes,
                 ha="left", va="center", fontsize=COLORBAR_LABEL_FONTSIZE)
    cbar.ax.xaxis.set_ticks_position("top")
    cbar.ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    cbar.ax.tick_params(axis="x", direction="in", labelsize=COLORBAR_TICK_FONTSIZE)

    draw_bottom_panel(ax_bottom, velocity_cbar_ax, t_edges, z_edges, delta_ebv_map,
                      cuts, clusters, ebv_vmin, ebv_vmax)
    # Mark the exact lower-panel field of view in the overview's data coordinates.
    zoom_xmin, zoom_xmax = TZ_BOTTOM_XLIM_PC
    zoom_zmin, zoom_zmax = TZ_BOTTOM_YLIM_PC
    zoom_box = Rectangle(
        (zoom_xmin, zoom_zmin), zoom_xmax - zoom_xmin, zoom_zmax - zoom_zmin,
        facecolor="none", edgecolor="0.25", linewidth=1.2,
        linestyle=(0, (5, 3)), zorder=19,
        path_effects=[pe.withStroke(linewidth=2.4, foreground="white")],
    )
    ax_top.add_patch(zoom_box)
    for overview_x, detail_x in ((zoom_xmin, 0.0), (zoom_xmax, 1.0)):
        # Draw above the overview background so the lines reach the box corners.
        ax_top.add_artist(ConnectionPatch(
            xyA=(overview_x, zoom_zmin), coordsA="data", axesA=ax_top,
            xyB=(detail_x, 1.0), coordsB="axes fraction", axesB=ax_bottom,
            color="0.5", linewidth=1.0, linestyle=(0, (5, 3)),
            clip_on=False, zorder=18,
        ))
    fig.savefig(OUT_FIG, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

def load_projected_hunt_konietzka_clusters(theta: float, point_t0: np.ndarray, intercept_kpc: float, half_width_pc: float):
    """Select complete young-cluster rows and project positions/LSR velocities.

    The ds normal offset and s/T coordinate are in pc; vs=vx*cos+vy*sin
    and vn=-vx*sin+vy*cos are in km/s. Keep ages strictly below 30 Myr,
    |ds| strictly below the band half-width, and positions inside the zoom.
    The raw Konietzka path constant is not read by this function."""
    df = pd.read_csv(HUNT_KONIETZKA_MINIMAL_CSV)
    needed = [
        "x_pc", "y_pc", "z_pc",
        "vx_lsr_kms", "vy_lsr_kms", "vz_lsr_kms",
        "age_myr", "source_key",
    ]
    df = df.dropna(subset=needed).copy()
    df = df[df["age_myr"] < BOTTOM_CLUSTER_AGE_CUT_MYR].copy()
    x = df["x_pc"].to_numpy()
    y = df["y_pc"].to_numpy()
    vx = df["vx_lsr_kms"].to_numpy()
    vy = df["vy_lsr_kms"].to_numpy()
    intercept_pc = 1000.0 * intercept_kpc

    ds = -np.sin(theta) * x + np.cos(theta) * y - intercept_pc * np.cos(theta)
    s_proj = (x / 1000.0 - point_t0[0]) * np.cos(theta) + (y / 1000.0 - point_t0[1]) * np.sin(theta)
    s_proj *= 1000.0
    vs = vx * np.cos(theta) + vy * np.sin(theta)
    vn = -vx * np.sin(theta) + vy * np.cos(theta)
    df = df.assign(
        x=x,
        y=y,
        z=df["z_pc"].to_numpy(),
        vx=vx,
        vy=vy,
        vz=df["vz_lsr_kms"].to_numpy(),
        age=df["age_myr"].to_numpy(),
        s=s_proj,
        ds=ds,
        vs=vs,
        vn=vn,
    )
    in_slice = np.abs(df["ds"]) < half_width_pc
    in_axes = df["s"].between(*TZ_BOTTOM_XLIM_PC) & df["z"].between(*TZ_BOTTOM_YLIM_PC)
    return df[in_slice & in_axes].copy()


def age_to_sizes(age: np.ndarray):
    """Map the selected age range linearly to marker areas from 1 to 40.

    Older clusters receive larger markers; equal ages receive the midpoint.
    This scaling is relative to the plotted sample, not a physical size."""
    s_min_size, s_max_size = 1, 40
    age_min, age_max = np.nanmin(age), np.nanmax(age)
    if age_max > age_min:
        return s_min_size + (s_max_size - s_min_size) * (age - age_min) / (age_max - age_min)
    return np.full_like(age, (s_min_size + s_max_size) / 2)


def draw_bottom_panel(ax, cbar_ax, t_edges, z_edges, delta_ebv_map, cuts, clusters, ebv_vmin: float, ebv_vmax: float):
    """Overlay young-cluster velocities on the grayscale dust zoom.

    Omit SB17 only from this panel's shell overlay. Marker shapes distinguish
    sources; marker sizes encode age, and both arrows and points use the
    same normal-velocity color scale. Arrow components lie in the T-Z plane
    and the quiver key supplies the 10 km/s display reference."""
    ax.pcolormesh(1000.0 * t_edges, 1000.0 * z_edges, delta_ebv_map, cmap="gray_r",
                  vmin=ebv_vmin, vmax=ebv_vmax, shading="auto", rasterized=True)
    ax.set_aspect("equal", adjustable="box")
    cuts_without_sb17 = cuts[~cuts["id"].eq(17)].copy()
    draw_superbubble_cuts_pc(ax, cuts_without_sb17, TZ_BOTTOM_XLIM_PC, TZ_BOTTOM_YLIM_PC, outline=True)
    draw_orion_eridanus_marker(ax)

    age = clusters["age"].to_numpy()
    sizes = age_to_sizes(age)
    for i, age_val in enumerate([5, 10, 20]):
        highlight_s = age_to_sizes(np.array([np.nanmin(age), age_val, np.nanmax(age)]))[1]
        ax.scatter(2000 * 0.95, 400 * (1 - 0.05 * i), s=highlight_s,
                   c="k", edgecolors="black", linewidth=1.2, label=f"Age = {age_val} Myr",
                   zorder=22)
    age_legend = ax.legend(loc="lower right", fontsize=8, labelspacing=0.35, handletextpad=0.5,
                           borderpad=0.35, framealpha=0.95)
    age_legend.set_zorder(30)

    cmap = mpl.colormaps.get_cmap("coolwarm")
    norm = Normalize(vmin=-20, vmax=20)
    sct = None
    source_markers = {
        "konietzka2023": "*",
        "hunt_stratified": "o",
    }
    for source_key, marker in source_markers.items():
        mask = clusters["source_key"].eq(source_key).to_numpy()
        if not np.any(mask):
            continue
        sct = ax.scatter(
            clusters.loc[mask, "s"].to_numpy(), clusters.loc[mask, "z"].to_numpy(),
            c=clusters.loc[mask, "vn"].to_numpy(), s=sizes[mask], cmap=cmap, norm=norm,
            marker=marker, alpha=0.6, zorder=20,
        )
    if sct is None:
        sct = ax.scatter(
            clusters["s"].to_numpy(), clusters["z"].to_numpy(),
            c=clusters["vn"].to_numpy(), s=sizes, cmap=cmap, norm=norm,
            alpha=0.6, zorder=20,
        )
    cbar_vn = ax.figure.colorbar(sct, cax=cbar_ax, orientation="horizontal")
    cbar_vn.ax.tick_params(axis="x", direction="in", labelsize=COLORBAR_TICK_FONTSIZE)
    cbar_vn.set_label(r"$v_n$ (km/s)", fontsize=COLORBAR_LABEL_FONTSIZE)

    q = ax.quiver(
        clusters["s"].to_numpy(), clusters["z"].to_numpy(),
        clusters["vs"].to_numpy(), clusters["vz"].to_numpy(), clusters["vn"].to_numpy(),
        cmap=cmap, norm=norm, angles="xy", scale_units="xy",
        scale=0.08, width=0.003, zorder=21,
    )
    ax.quiverkey(q, X=0.1, Y=0.9, U=10, label="10 km/s", labelpos="E", coordinates="axes")

    ax.set_xlim(*TZ_BOTTOM_XLIM_PC)
    ax.set_ylim(*TZ_BOTTOM_YLIM_PC)
    ax.text(0.01, 0.94, "(b)", transform=ax.transAxes, ha="left", va="top",
            fontsize=PANEL_LABEL_FONTSIZE, fontweight="bold",
            path_effects=[pe.withStroke(linewidth=3.0, foreground="white")])
    ax.set_xlabel("T (pc)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel("Z (pc)", fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(axis="both", which="major", labelsize=TICK_LABEL_FONTSIZE)
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True)
    ax.grid(True, linestyle="--", alpha=0.5)


def main():
    """Select/cache the dust band, compute central-plane shell cuts, and render both panels."""
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FINAL_FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Computing RW cross-sections at angle={args.angle:g} deg, intercept={args.intercept:g} kpc.")
    theta, unit_t, unit_n, point_t0 = line_geometry(args.angle, args.intercept)
    default_rw_slice = (
        abs(args.angle - ANGLE_DEG) < 1e-9
        and abs(args.intercept - INTERCEPT_KPC) < 1e-9
        and abs(args.hw - HALF_WIDTH_KPC) < 1e-9
    )
    if default_rw_slice:
        slice_df = load_or_create_rw_dust_slice()
        print(f"Loaded RW dust slice with {len(slice_df)} cells.")
    else:
        dust_df = load_dust_subset()
        print(f"Loaded 3D dust subset with {len(dust_df)} cells.")
        slice_df = build_band_slice_df(dust_df, theta, point_t0, args.intercept, args.hw)

    print("Computing OSB intersections with the RW-aligned slice.")
    superbubbles = load_superbubbles()
    # Shell outlines use the central plane, even though dust/clusters use a band.
    # In particular, SB31 follows the original fit without shrinking its boundary.
    cuts = compute_all_cuts(superbubbles, unit_t, point_t0, unit_n=unit_n, half_width=0.0)
    print(f"Computed projected cuts for {len(cuts)} OSBs.")

    print("Building top-panel dust map.")
    t_edges, z_edges, mean_dust = build_dust_tz_map_from_slice(slice_df)
    slice_thickness_kpc = 2.0 * args.hw
    # Mean differential reddening [mag/kpc] times slab thickness [kpc].
    # This assumes the sampled mean represents the full depth of each T-Z bin.
    delta_ebv_map = mean_dust * slice_thickness_kpc
    ebv_vmin = DELTA_EBV_VMIN
    ebv_vmax = DELTA_EBV_VMAX

    print("Loading projected young clusters for the bottom panel.")
    clusters = load_projected_hunt_konietzka_clusters(theta, point_t0, args.intercept, args.hw * 1000.0)
    source_counts = clusters["source_key"].value_counts().to_dict()
    print(
        f"Projected young-cluster sample count: {len(clusters)}; source: {source_counts}"
    )
    print("Rendering combined RW dust and young-cluster figure.")
    draw_figure(t_edges, z_edges, delta_ebv_map, cuts, clusters, args.dpi, ebv_vmin, ebv_vmax)
    print(f"  saved: {OUT_FIG}")

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
