#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Script 0a: 3D dust-map data preparation

Purpose
-------
Build raw/smoothed 3D dust cubes, XY mean dust map, and Radcliffe-Wave dust slice from dustmaps3d.

Method overview
---------------
1. Build the Cartesian dust grid and convert it to Galactic coordinates.
2. Query dustmaps3d, mask invalid voxels, and build raw and smoothed dust products.
3. Export shared dust-map products used by later geometry, plotting, and visualization
   stages.

Main inputs
-----------
- dustmaps3d local installation and Wang et al. (2025) dust-map cache

Main outputs
------------
- Shared intermediate inputs in ../results/intermediate_output/3d_dust_map_products/ for scripts 5, 8, 11, 12, 23, 24, 25.

Figure/table role
-----------------
Shared intermediate inputs in ../results/intermediate_output/3d_dust_map_products/ for scripts 5, 8, 11, 12, 23, 24, 25.

Runtime and data notes
----------------------
Reference runtime: 20 min 58.49 s. The current implementation uses X-slab queries, memmap temporary cubes, and chunked Parquet writes so it can run inside capsule environments with tighter RAM limits. First use of dustmaps3d may also download the 3D dust-map cache. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

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

from pathlib import Path
import tempfile
from typing import Tuple
import time

import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord, CartesianRepresentation
from dustmaps3d import dustmaps3d
from scipy.ndimage import gaussian_filter, distance_transform_edt
import pyarrow as pa
import pyarrow.parquet as pq


# ------------------------------
# Default parameters
# ------------------------------
GRID_X = (-4.5, 4.5)
GRID_Y = (-4.5, 4.5)
GRID_Z = (-1.0, 1.0)
STEP = 0.01  # kpc

SMOOTH_SIGMA = 2  # Gaussian sigma (grid cells)
QUERY_X_BLOCK = 8  # X-slab width for memory-bounded dust-map queries.
SMOOTH_X_BLOCK = 16  # X-slab width for memory-bounded smoothing and export.

# Radcliffe Wave reference geometry and slice configuration.
BAND_ANGLE_DEG = 60.0
BAND_OFFSET_KPC = 0.6
BAND_HALF_WIDTH_KPC = 0.1

OUT_DIR = Path("..") / "results" / "intermediate_output" / "3d_dust_map_products"
OUT_RAW = OUT_DIR / "raw_3d_dust_cube.parquet"
OUT_SMOOTH = OUT_DIR / "smoothed_3d_dust_cube.parquet"
OUT_2D = OUT_DIR / "xy_mean_dust_map.parquet"
OUT_S = OUT_DIR / "radcliffe_wave_dust_slice.csv"


def make_axis(bounds: Tuple[float, float], step: float) -> np.ndarray:
    """Return cell-centred grid coordinates for one Cartesian axis."""
    n_cells = int(round((bounds[1] - bounds[0]) / step))
    return (bounds[0] + step / 2 + np.arange(n_cells, dtype=np.float32) * step).astype(np.float32)


def flattened_coordinate_block(
    xv: np.ndarray,
    y_base: np.ndarray,
    z_base: np.ndarray,
    yz_size: int,
    x_start: int,
    x_stop: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return flattened coordinates for one X slab in original meshgrid order."""
    x_chunk = xv[x_start:x_stop]
    n_x = len(x_chunk)
    x = np.repeat(x_chunk, yz_size).astype(np.float32, copy=False)
    y = np.tile(y_base, n_x).astype(np.float32, copy=False)
    z = np.tile(z_base, n_x).astype(np.float32, copy=False)
    return x, y, z


def to_galactic(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert Sun-centered Cartesian (kpc) to Galactic (l,b,d)."""
    rep = CartesianRepresentation(x * u.kpc, y * u.kpc, z * u.kpc)
    sky = SkyCoord(rep, frame="galactic")
    return sky.l.deg.astype(np.float32), sky.b.deg.astype(np.float32), sky.distance.kpc.astype(np.float32)


def query_dust(l: np.ndarray, b: np.ndarray, d: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Query dustmaps3d: returns (dust_raw, d_max) arrays."""
    _, dust_raw, _, dmax = dustmaps3d(l, b, d)
    return np.asarray(dust_raw, dtype=np.float32), np.asarray(dmax, dtype=np.float32)


def write_float32_parquet_chunk(
    writer: pq.ParquetWriter | None,
    path: Path,
    names: Tuple[str, ...],
    columns: Tuple[np.ndarray, ...],
) -> pq.ParquetWriter:
    """Append one float32 chunk to a Parquet file and return the active writer."""
    table = pa.Table.from_arrays(
        [pa.array(np.asarray(col, dtype=np.float32), type=pa.float32()) for col in columns],
        names=list(names),
    )
    if writer is None:
        writer = pq.ParquetWriter(path, table.schema, compression="snappy")
    writer.write_table(table)
    return writer


def query_and_write_raw_cube(
    xv: np.ndarray,
    yv: np.ndarray,
    zv: np.ndarray,
    raw_mm: np.memmap,
) -> None:
    """Query dustmaps3d in X slabs, write raw Parquet, and populate raw memmap."""
    y_base = np.repeat(yv, len(zv)).astype(np.float32, copy=False)
    z_base = np.tile(zv, len(yv)).astype(np.float32, copy=False)
    yz_size = len(y_base)
    writer = None
    try:
        for x_start in range(0, len(xv), QUERY_X_BLOCK):
            x_stop = min(len(xv), x_start + QUERY_X_BLOCK)
            x, y, z = flattened_coordinate_block(xv, y_base, z_base, yz_size, x_start, x_stop)
            l, b, d = to_galactic(x, y, z)
            dust_raw, dmax = query_dust(l, b, d)
            dust_masked = np.where(d <= dmax, dust_raw, np.nan).astype(np.float32)
            finite_mask = ~np.isnan(dust_masked)
            dust_masked[(finite_mask) & (dust_masked < 1e-4)] = 0.0
            raw_mm[x_start:x_stop, :, :] = dust_masked.reshape((x_stop - x_start, len(yv), len(zv)))
            writer = write_float32_parquet_chunk(
                writer,
                OUT_RAW,
                ("X", "Y", "Z", "dust_raw"),
                (x, y, z, dust_masked),
            )
            print(f"Queried and wrote raw dust X slabs {x_start}:{x_stop} / {len(xv)}", flush=True)
    finally:
        if writer is not None:
            writer.close()
    raw_mm.flush()


def build_nearest_filled_memmap(raw_mm: np.memmap, nan_mask: np.ndarray, tmpdir: Path) -> np.memmap:
    """Create a nearest-neighbor-filled dust cube without keeping EDT indices in RAM."""
    if not np.any(nan_mask):
        return raw_mm

    nx, ny, nz = raw_mm.shape
    idx_path = tmpdir / "nearest_indices_i4.dat"
    filled_path = tmpdir / "nearest_filled_f4.dat"
    idx = np.memmap(idx_path, dtype=np.int32, mode="w+", shape=(3, nx, ny, nz))
    distance_transform_edt(
        nan_mask,
        return_distances=False,
        return_indices=True,
        indices=idx,
    )
    idx.flush()

    filled_mm = np.memmap(filled_path, dtype=np.float32, mode="w+", shape=raw_mm.shape)
    for x_start in range(0, nx, SMOOTH_X_BLOCK):
        x_stop = min(nx, x_start + SMOOTH_X_BLOCK)
        filled_mm[x_start:x_stop, :, :] = raw_mm[
            idx[0, x_start:x_stop, :, :],
            idx[1, x_start:x_stop, :, :],
            idx[2, x_start:x_stop, :, :],
        ]
        print(f"Nearest-filled X slabs {x_start}:{x_stop} / {nx}", flush=True)
    filled_mm.flush()
    del idx
    return filled_mm


def smooth_cube_blockwise(
    filled_mm: np.memmap,
    nan_mask: np.ndarray,
    smooth_mm: np.memmap,
) -> None:
    """Apply the original nearest-fill plus Gaussian smoothing in X slabs."""
    nx = filled_mm.shape[0]
    radius = int(np.ceil(4.0 * float(SMOOTH_SIGMA)))
    for x_start in range(0, nx, SMOOTH_X_BLOCK):
        x_stop = min(nx, x_start + SMOOTH_X_BLOCK)
        halo_start = max(0, x_start - radius)
        halo_stop = min(nx, x_stop + radius)
        block = np.asarray(filled_mm[halo_start:halo_stop, :, :], dtype=np.float32)
        smooth_block = gaussian_filter(block, sigma=SMOOTH_SIGMA, mode="reflect", truncate=4.0)
        inner_start = x_start - halo_start
        inner_stop = inner_start + (x_stop - x_start)
        out_block = smooth_block[inner_start:inner_stop, :, :].astype(np.float32, copy=True)
        out_block[nan_mask[x_start:x_stop, :, :]] = np.nan
        smooth_mm[x_start:x_stop, :, :] = out_block
        print(f"Smoothed X slabs {x_start}:{x_stop} / {nx}", flush=True)
    smooth_mm.flush()


def write_smoothed_cube_parquet(
    xv: np.ndarray,
    yv: np.ndarray,
    zv: np.ndarray,
    smooth_mm: np.memmap,
) -> None:
    """Write smoothed dust memmap to Parquet in original meshgrid row order."""
    y_base = np.repeat(yv, len(zv)).astype(np.float32, copy=False)
    z_base = np.tile(zv, len(yv)).astype(np.float32, copy=False)
    yz_size = len(y_base)
    writer = None
    try:
        for x_start in range(0, len(xv), SMOOTH_X_BLOCK):
            x_stop = min(len(xv), x_start + SMOOTH_X_BLOCK)
            x, y, z = flattened_coordinate_block(xv, y_base, z_base, yz_size, x_start, x_stop)
            dust = np.asarray(smooth_mm[x_start:x_stop, :, :], dtype=np.float32).ravel()
            writer = write_float32_parquet_chunk(
                writer,
                OUT_SMOOTH,
                ("X", "Y", "Z", "dust"),
                (x, y, z, dust),
            )
            print(f"Wrote smoothed dust X slabs {x_start}:{x_stop} / {len(xv)}", flush=True)
    finally:
        if writer is not None:
            writer.close()


def write_xy_mean_map(xv: np.ndarray, yv: np.ndarray, zv: np.ndarray, smooth_mm: np.memmap) -> None:
    """Write the |Z| <= 0.3 kpc XY mean map without materializing the full cube."""
    z_mask = np.abs(zv) <= 0.3
    plane = np.asarray(smooth_mm[:, :, z_mask], dtype=np.float32)
    valid = np.isfinite(plane)
    count = valid.sum(axis=2)
    sum_map = np.where(valid, plane, 0.0).sum(axis=2, dtype=np.float32)
    mean_map = np.full(count.shape, np.nan, dtype=np.float32)
    np.divide(sum_map, count, out=mean_map, where=count > 0)
    X2, Y2 = np.meshgrid(xv, yv, indexing="ij")
    table = pa.Table.from_arrays(
        [
            pa.array(X2.ravel().astype(np.float32), type=pa.float32()),
            pa.array(Y2.ravel().astype(np.float32), type=pa.float32()),
            pa.array(mean_map.ravel().astype(np.float32), type=pa.float32()),
        ],
        names=["X", "Y", "dust"],
    )
    pq.write_table(table, OUT_2D, compression="snappy")


def write_radcliffe_wave_slice(xv: np.ndarray, yv: np.ndarray, zv: np.ndarray, smooth_mm: np.memmap) -> None:
    """Write the Radcliffe-Wave band-limited dust slice in chunks."""
    theta = np.deg2rad(BAND_ANGLE_DEG)
    x0 = -BAND_OFFSET_KPC * np.sin(theta) * np.cos(theta)
    y0 = BAND_OFFSET_KPC * (np.cos(theta) ** 2)
    y_base = np.repeat(yv, len(zv)).astype(np.float32, copy=False)
    z_base = np.tile(zv, len(yv)).astype(np.float32, copy=False)
    yz_size = len(y_base)

    wrote_header = False
    for x_start in range(0, len(xv), SMOOTH_X_BLOCK):
        x_stop = min(len(xv), x_start + SMOOTH_X_BLOCK)
        x, y, z = flattened_coordinate_block(xv, y_base, z_base, yz_size, x_start, x_stop)
        distances = -np.sin(theta) * x + np.cos(theta) * y - BAND_OFFSET_KPC * np.cos(theta)
        mask_band = np.abs(distances) <= BAND_HALF_WIDTH_KPC
        if not np.any(mask_band):
            continue
        dust = np.asarray(smooth_mm[x_start:x_stop, :, :], dtype=np.float32).ravel()
        s = ((x - x0) * np.cos(theta) + (y - y0) * np.sin(theta)).astype(np.float32)
        df_band = pd.DataFrame(
            {
                "X": x[mask_band],
                "Y": y[mask_band],
                "Z": z[mask_band],
                "dust": dust[mask_band],
                "s": s[mask_band],
            }
        )
        df_band.to_csv(
            OUT_S,
            mode="w" if not wrote_header else "a",
            header=not wrote_header,
            index=False,
            float_format="%.8g",
        )
        wrote_header = True
    if not wrote_header:
        pd.DataFrame(columns=["X", "Y", "Z", "dust", "s"]).to_csv(OUT_S, index=False)


def main() -> None:
    """Run Script 0a from validated inputs to the documented outputs."""
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    xv = make_axis(GRID_X, STEP)
    yv = make_axis(GRID_Y, STEP)
    zv = make_axis(GRID_Z, STEP)
    nx, ny, nz = len(xv), len(yv), len(zv)
    print(f"Preparing dust grid: {nx} x {ny} x {nz} cells", flush=True)

    # The science product is unchanged from the original implementation, but
    # Code Ocean capsules can have tighter memory limits than local workstations.
    # Temporary memmaps and chunked Parquet row groups avoid materializing all
    # 1.62e8 coordinates and dust values as in-memory DataFrames.
    with tempfile.TemporaryDirectory(prefix="osb_dust_") as tmp:
        tmpdir = Path(tmp)
        raw_mm = np.memmap(tmpdir / "raw_dust_f4.dat", dtype=np.float32, mode="w+", shape=(nx, ny, nz))
        smooth_mm = np.memmap(tmpdir / "smooth_dust_f4.dat", dtype=np.float32, mode="w+", shape=(nx, ny, nz))

        # 1) Query dust and save the raw 3D cube.
        query_and_write_raw_cube(xv, yv, zv, raw_mm)
        print(f"Saved raw 3D cube: {OUT_RAW}")

        # 2) Smooth 3D cube with the original nearest-fill convention.
        nan_mask = np.isnan(raw_mm)
        filled_mm = build_nearest_filled_memmap(raw_mm, nan_mask, tmpdir)
        smooth_cube_blockwise(filled_mm, nan_mask, smooth_mm)

        # 3) Save smoothed 3D cube and derived 2D/slice products.
        write_smoothed_cube_parquet(xv, yv, zv, smooth_mm)
        print(f"Saved smoothed 3D cube: {OUT_SMOOTH}")

        write_xy_mean_map(xv, yv, zv, smooth_mm)
        print(f"Saved 2D mean map: {OUT_2D}")

        write_radcliffe_wave_slice(xv, yv, zv, smooth_mm)
        print(f"Saved RW along-path dust slice: {OUT_S}")

        del filled_mm
        del smooth_mm
        del raw_mm

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
