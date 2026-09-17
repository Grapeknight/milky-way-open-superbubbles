#!/usr/bin/env python
"""Script 0b: CO and HI PV-cache preparation

Purpose
-------
Export compact MWISP/CfA/HI4PI PV caches from external original FITS cubes. This step cannot be run from the package alone because the several-hundred-GB raw cubes are not included.

Method overview
---------------
1. Validate that the external survey FITS cubes are available to the local user.
2. Derive longitude, latitude, and velocity axes from FITS headers and collapse the data
   into compact PV products.
3. Write only the small cache files consumed by the downstream PV plotting scripts.

Main inputs
-----------
- external MWISP 12CO FITS cube (not included)
- external CfA CO FITS cube (not included)
- external HI4PI HI FITS cube (not included)

Main outputs
------------
- ../data/COdata/, ../data/HIdata/; these compact caches are included and used by scripts 21 and 22.

Figure/table role
-----------------
../data/COdata/, ../data/HIdata/; these compact caches are included and used by scripts 21 and 22.

Runtime and data notes
----------------------
External raw-data export timing estimate in README.md: > 6 h. This script was not included in the standard reproducible run because the several-hundred-GB raw survey cubes are not packaged. All normal paths are resolved relative to the code directory using the current capsule tree: ../data for inputs and ../results for generated outputs unless an external-data caveat below says otherwise.

External-data caveat
--------------------
The several-hundred-GB MWISP, CfA, and HI4PI raw cubes are not included; this file documents how the compact PV caches were exported from those raw products.

Reading the code
----------------
Start with the path and scientific settings below, then follow main() at
the end of the file. The method overview above describes the order of the
analysis steps; the input/output lists identify the upstream data and
products written by this stage. Relative data paths are resolved from code/.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path
from typing import Iterable

import numpy as np
from astropy.io import fits

try:
    from tqdm import tqdm
except ImportError:  # pragma: no cover - convenience fallback for lean envs.
    tqdm = None


RAW_SURVEY_DIR = Path(os.environ.get("OSB_RAW_SURVEY_DIR", "external_raw_data"))
DEFAULT_MWISP_FITS = Path(
    os.environ.get("MWISP_DATA_PATH", str(RAW_SURVEY_DIR / "MWISP" / "mosaic_12CO.fits"))
)
DEFAULT_CFA_FITS = Path(
    os.environ.get("CFA_CO_DATA_PATH", str(RAW_SURVEY_DIR / "CfA_CO" / "DT22_DHT_masked.fits"))
)
DEFAULT_HI4PI_FITS = Path(
    os.environ.get("HI4PI_DATA_PATH", str(RAW_SURVEY_DIR / "HI4PI" / "CAR.fits"))
)

MWISP_OUT_NAME = "MWISP_12CO_PV_cache.npz"
CFA_DEFAULT_OUT_NAME = "CfA_12CO_PV_cache.npz"
CFA_OUT_PREFIX = "CfA_12CO_PV_babs"
HI_DEFAULT_OUT_NAME = "HI4PI_HI_PV_cache.npy"
HI_OUT_PREFIX = "HI4PI_HI_PV_babs"
HI_AXES_NAME = "HI4PI_HI_PV_axes.npz"
HI_SUMMARY_NAME = "HI4PI_HI_PV_cache_summary.json"


def default_base_dir() -> Path:
    script_dir = Path(__file__).resolve().parent
    return script_dir.parent if script_dir.name.lower() == "code" else script_dir


def parse_args() -> argparse.Namespace:
    """Build the command-line interface for Script 0b."""
    parser = argparse.ArgumentParser(
        description="Build CO and HI PV caches under ../data/COdata and ../data/HIdata."
    )
    parser.add_argument(
        "--products",
        nargs="+",
        default=["all"],
        choices=["all", "co", "hi", "mwisp", "cfa", "hi4pi"],
        help="Products to build. Default: all.",
    )
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=default_base_dir(),
        help="Project root. Default: the directory containing this script.",
    )
    parser.add_argument("--co-out-dir", type=Path, default=None, help="CO cache output directory.")
    parser.add_argument("--hi-out-dir", type=Path, default=None, help="HI cache output directory.")
    parser.add_argument("--mwisp-fits", type=Path, default=DEFAULT_MWISP_FITS, help="MWISP mosaic_12CO.fits.")
    parser.add_argument("--cfa-fits", type=Path, default=DEFAULT_CFA_FITS, help="CfA DT22+DHT FITS cube.")
    parser.add_argument("--hi4pi-fits", type=Path, default=DEFAULT_HI4PI_FITS, help="HI4PI CAR.fits cube.")
    parser.add_argument(
        "--cfa-b-limits",
        nargs="+",
        default=["15"],
        help="CfA latitude cuts in degrees. Accepts space-separated or comma-separated values. Default: 15.",
    )
    parser.add_argument(
        "--hi-b-max-deg",
        type=float,
        default=15.0,
        help="HI4PI absolute Galactic latitude cut for PV integration. Default: 15.",
    )
    parser.add_argument(
        "--mwisp-lat-chunk",
        type=int,
        default=1,
        help="MWISP latitude rows processed per chunk. Default: 1.",
    )
    parser.add_argument(
        "--hi4pi-vel-chunk",
        type=int,
        default=8,
        help="HI4PI velocity channels processed per chunk. Default: 8.",
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Do not rebuild products whose output cache already exists.",
    )
    return parser.parse_args()


def axis_from_header(header: fits.Header, axis: int) -> np.ndarray:
    n = int(header[f"NAXIS{axis}"])
    crpix = float(header[f"CRPIX{axis}"])
    crval = float(header[f"CRVAL{axis}"])
    cdelt = float(header[f"CDELT{axis}"])
    return crval + (np.arange(n, dtype=np.float64) + 1.0 - crpix) * cdelt


def progress(iterable: Iterable[int], **kwargs: object) -> Iterable[int]:
    if tqdm is None:
        return iterable
    return tqdm(iterable, **kwargs)


def require_existing_file(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} not found: {path}")
    if not path.is_file():
        raise FileNotFoundError(f"{label} is not a file: {path}")


def elapsed_text(seconds: float) -> str:
    minutes, sec = divmod(int(max(seconds, 0.0)), 60)
    hours, minutes = divmod(minutes, 60)
    return f"{hours:02d}:{minutes:02d}:{sec:02d}"


def babs_suffix(b_limit: float) -> str:
    return str(int(b_limit)) if float(b_limit).is_integer() else str(b_limit).replace(".", "p")


def hi_b_tag(b_max_deg: float) -> str:
    return f"{b_max_deg:04.1f}".replace(".", "p")


def parse_float_list(values: list[str]) -> list[float]:
    parsed: list[float] = []
    for value in values:
        parsed.extend(float(part) for part in value.split(",") if part.strip())
    if not parsed:
        raise ValueError("At least one latitude cut is required.")
    return parsed


def selected_products(raw_products: list[str]) -> set[str]:
    selected: set[str] = set()
    for product in raw_products:
        if product == "all":
            selected.update({"mwisp", "cfa", "hi4pi"})
        elif product == "co":
            selected.update({"mwisp", "cfa"})
        elif product == "hi":
            selected.add("hi4pi")
        else:
            selected.add(product)
    return selected


def build_mwisp_lv_cache(
    fits_path: Path,
    out_dir: Path,
    lat_chunk: int = 1,
    skip_existing: bool = False,
) -> Path:
    out_path = out_dir / MWISP_OUT_NAME
    if skip_existing and out_path.exists():
        print(f"Skip existing MWISP cache: {out_path}", flush=True)
        return out_path

    require_existing_file(fits_path, "MWISP FITS")
    out_dir.mkdir(parents=True, exist_ok=True)
    lat_chunk = max(int(lat_chunk), 1)

    with fits.open(fits_path, mode="readonly", memmap=True) as hdul:
        header = hdul[0].header
        data = hdul[0].data

        lon = axis_from_header(header, 1)
        lat = axis_from_header(header, 2)
        vel = axis_from_header(header, 3) / 1000.0
        db_deg = abs(float(header["CDELT2"]))
        bmin = float(lat.min())
        bmax = float(lat.max())

        nv = int(header["NAXIS3"])
        nx = int(header["NAXIS1"])
        ny = int(header["NAXIS2"])
        lv_map = np.zeros((nv, nx), dtype=np.float32)

        starts = range(0, ny, lat_chunk)
        total = (ny + lat_chunk - 1) // lat_chunk
        start_time = time.perf_counter()
        print(f"Building MWISP 12CO PV cache from {fits_path}", flush=True)

        for y0 in progress(starts, total=total, desc="MWISP latitude chunks"):
            y1 = min(y0 + lat_chunk, ny)
            block = np.asarray(data[:, y0:y1, :], dtype=np.float32)
            lv_chunk = np.nansum(block, axis=1, dtype=np.float64) * db_deg
            lv_map += lv_chunk.astype(np.float32)

        print(f"MWISP integration finished in {elapsed_text(time.perf_counter() - start_time)}", flush=True)

    lv_map[lv_map < 0] = 0
    np.savez_compressed(
        out_path,
        lv_map=lv_map.astype(np.float32),
        lon=lon.astype(np.float64),
        vel=vel.astype(np.float64),
        bmin=np.array([bmin], dtype=np.float32),
        bmax=np.array([bmax], dtype=np.float32),
    )
    print(f"Saved MWISP 12CO PV cache: {out_path}", flush=True)
    return out_path


def build_cfa_lv_cache(
    fits_path: Path,
    out_dir: Path,
    b_limit: float,
    skip_existing: bool = False,
) -> Path:
    suffix = babs_suffix(b_limit)
    out_name = CFA_DEFAULT_OUT_NAME if float(b_limit) == 15.0 else f"{CFA_OUT_PREFIX}{suffix}_cache.npz"
    out_path = out_dir / out_name
    if skip_existing and out_path.exists():
        print(f"Skip existing CfA 12CO PV cache: {out_path}", flush=True)
        return out_path

    require_existing_file(fits_path, "CfA FITS")
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"Building CfA 12CO PV cache |b| <= {b_limit:g} deg from {fits_path}", flush=True)
    with fits.open(fits_path, memmap=True) as hdul:
        header = hdul[0].header
        cube = hdul[0].data
        velocity = axis_from_header(header, 1)
        longitude = np.mod(axis_from_header(header, 2), 360.0)
        latitude = axis_from_header(header, 3)

        order = np.argsort(longitude)
        lon_sorted = longitude[order]
        lat_mask = np.abs(latitude) <= b_limit
        if not np.any(lat_mask):
            raise ValueError(f"No CfA latitude pixels found within |b| <= {b_limit:g} deg.")

        dlat = abs(float(header["CDELT3"]))
        lv_map = np.nansum(cube[lat_mask, :, :][:, order, :], axis=0, dtype=np.float64) * dlat
        lv_map = lv_map.T.astype(np.float32)
        lv_map[lv_map < 0] = 0

    np.savez_compressed(
        out_path,
        lv_map=lv_map,
        lon=lon_sorted.astype(np.float64),
        vel=velocity.astype(np.float64),
        b_limit=np.array([b_limit], dtype=float),
    )
    print(f"Saved CfA 12CO PV cache: {out_path}", flush=True)
    return out_path


def write_hi_axes_and_summary(
    out_dir: Path,
    source_path: Path,
    header: fits.Header,
    lon_deg: np.ndarray,
    lat_deg: np.ndarray,
    vel_kms: np.ndarray,
    b_max_deg: float,
    lv_path: Path,
) -> None:
    np.savez_compressed(
        out_dir / HI_AXES_NAME,
        lon_deg=lon_deg.astype(np.float32),
        lat_deg=lat_deg.astype(np.float32),
        vel_kms=vel_kms.astype(np.float32),
    )

    summary = {
        "package_purpose": "HI4PI HI PV cache package without the original CAR.fits cube",
        "source_fits": str(source_path),
        "packaged_build_script": Path(__file__).name,
        "pv_cache_file": lv_path.name,
        "axes_file": HI_AXES_NAME,
        "latitude_cut_deg": float(b_max_deg),
        "shape_v_b_l": [
            int(header["NAXIS3"]),
            int(header["NAXIS2"]),
            int(header["NAXIS1"]),
        ],
        "pv_cache_shape": [
            int(header["NAXIS3"]),
            int(header["NAXIS1"]),
        ],
        "lon_axis_length": int(lon_deg.size),
        "lat_axis_length": int(lat_deg.size),
        "vel_axis_length": int(vel_kms.size),
        "ctype": [header.get("CTYPE1"), header.get("CTYPE2"), header.get("CTYPE3")],
        "cunit": [header.get("CUNIT1"), header.get("CUNIT2"), header.get("CUNIT3")],
        "bunit": header.get("BUNIT"),
        "lon_range_deg": [float(lon_deg.min()), float(lon_deg.max())],
        "lat_range_deg": [float(lat_deg.min()), float(lat_deg.max())],
        "vel_range_kms": [float(vel_kms.min()), float(vel_kms.max())],
        "dl_deg": float(abs(header["CDELT1"])),
        "db_deg": float(abs(header["CDELT2"])),
        "dv_kms": float(abs(header["CDELT3"]) / 1000.0),
        "note": (
            "The PV cache is built by summing HI brightness temperature over "
            "latitude within |b|<=b_max_deg and multiplying by db in degrees."
        ),
    }
    (out_dir / HI_SUMMARY_NAME).write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )


def build_hi4pi_lv_cache(
    fits_path: Path,
    out_dir: Path,
    b_max_deg: float = 15.0,
    vel_chunk: int = 8,
    skip_existing: bool = False,
) -> Path:
    tag = hi_b_tag(b_max_deg)
    out_name = HI_DEFAULT_OUT_NAME if float(b_max_deg) == 15.0 else f"{HI_OUT_PREFIX}{tag}deg_cache.npy"
    out_path = out_dir / out_name
    if skip_existing and out_path.exists() and (out_dir / HI_AXES_NAME).exists():
        print(f"Skip existing HI4PI HI PV cache: {out_path}", flush=True)
        return out_path

    require_existing_file(fits_path, "HI4PI FITS")
    out_dir.mkdir(parents=True, exist_ok=True)
    vel_chunk = max(int(vel_chunk), 1)

    with fits.open(fits_path, memmap=True) as hdul:
        header = hdul[0].header.copy()
        data = hdul[0].data

        lon_deg = axis_from_header(header, 1)
        lat_deg = axis_from_header(header, 2)
        vel_kms = axis_from_header(header, 3) / 1000.0
        write_hi_axes_and_summary(out_dir, fits_path, header, lon_deg, lat_deg, vel_kms, b_max_deg, out_path)

        db_deg = abs(float(header["CDELT2"]))
        nv = int(header["NAXIS3"])
        nx = int(header["NAXIS1"])
        lat_mask = np.abs(lat_deg) <= b_max_deg
        if not np.any(lat_mask):
            raise ValueError(f"No HI4PI latitude pixels found within |b| <= {b_max_deg:g} deg.")

        lv_map = np.zeros((nv, nx), dtype=np.float32)
        starts = range(0, nv, vel_chunk)
        total = (nv + vel_chunk - 1) // vel_chunk
        print(f"Building HI4PI HI PV cache |b| <= {b_max_deg:g} deg from {fits_path}", flush=True)
        for start in progress(starts, total=total, desc="HI4PI velocity chunks"):
            stop = min(start + vel_chunk, nv)
            chunk = data[start:stop, :, :]
            lv_chunk = np.nansum(chunk[:, lat_mask, :], axis=1, dtype=np.float64) * db_deg
            lv_map[start:stop, :] = lv_chunk.astype(np.float32)

    np.save(out_path, lv_map)
    print(f"Saved HI4PI HI PV cache: {out_path}", flush=True)
    print(f"Saved HI4PI axes: {out_dir / HI_AXES_NAME}", flush=True)
    print(f"Saved HI4PI summary: {out_dir / HI_SUMMARY_NAME}", flush=True)
    return out_path


def main() -> None:
    """Run Script 0b from validated inputs to the documented outputs."""
    args = parse_args()
    products = selected_products(args.products)
    base_dir = args.base_dir
    co_out_dir = args.co_out_dir or base_dir / "data" / "COdata"
    hi_out_dir = args.hi_out_dir or base_dir / "data" / "HIdata"

    if "mwisp" in products:
        build_mwisp_lv_cache(
            args.mwisp_fits,
            co_out_dir,
            lat_chunk=args.mwisp_lat_chunk,
            skip_existing=args.skip_existing,
        )

    if "cfa" in products:
        for b_limit in parse_float_list(args.cfa_b_limits):
            build_cfa_lv_cache(
                args.cfa_fits,
                co_out_dir,
                b_limit=b_limit,
                skip_existing=args.skip_existing,
            )

    if "hi4pi" in products:
        build_hi4pi_lv_cache(
            args.hi4pi_fits,
            hi_out_dir,
            b_max_deg=args.hi_b_max_deg,
            vel_chunk=args.hi4pi_vel_chunk,
            skip_existing=args.skip_existing,
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
