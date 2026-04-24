#!/usr/bin/env python3
"""
Register H&E and PSR stained mouse liver WSIs.

Sections are 2-3 apart, so a 3-stage Elastix pipeline is used:
  1. rigid  — correct slide placement (rotation + translation)
  2. affine — correct shear / global scale from cutting angle differences
  3. nl     — correct local tissue deformation from the section gap

Similarity metric: Mutual Information (multimodal-safe, works across H&E ↔ PSR).
Fixed image : PSR   (reference space)
Moving image: H&E   (warped into PSR space)

Usage:
    python register.py --he path/to/he.svs --psr path/to/psr.svs --out ./output --qc

    # If your scanner pixel size differs from the default 0.23 µm (40× Leica):
    python register.py --he he.svs --psr psr.svs --pixel-size 0.46  # 20× objective
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

try:
    from wsireg.wsireg2d import WsiReg2D
except ImportError:
    sys.exit(
        "wsireg not found.\n"
        "Install with: pip install wsireg\n"
        "If on Apple Silicon: conda install -c conda-forge itk-elastix first."
    )

try:
    import openslide
except ImportError:
    sys.exit(
        "openslide-python not found.\n"
        "Install with: pip install openslide-python\n"
        "You also need the OpenSlide C library: brew install openslide  (macOS)"
    )


# ─── Scanner / imaging defaults ───────────────────────────────────────────────

DEFAULT_PIXEL_SIZE_UM = 0.23   # µm/px at full resolution — 40× Leica Aperio
QC_PYRAMID_LEVEL      = 5      # which WSI pyramid level to use for QC thumbnails


# ─── Per-stain preprocessing passed to wsireg ────────────────────────────────
#
# Green channel (index 1) is used for both stains because:
#   H&E : nuclei and cytoplasm absorb green light → dark on bright background,
#          giving high-contrast tissue structure.
#   PSR : red collagen absorbs green → appears dark; hematoxylin counterstain
#          also shows nuclei dark on the green channel.
# Using the same channel type on both sides makes the mutual-information
# metric easier to optimise without losing the multimodal safety it provides.
#
# CLAHE (contrast-limited adaptive histogram equalisation) equalises local
# contrast so that pale pericentral zones and dense portal regions contribute
# equally to the registration — important for mouse liver lobular heterogeneity.

HE_PREPRO = {
    "as_uint8":         True,
    "ch_indices":       [1],      # green channel of RGB
    "contrast_enhance": "clahe",
}

PSR_PREPRO = {
    "as_uint8":         True,
    "ch_indices":       [1],      # green channel of RGB
    "contrast_enhance": "clahe",
}


# ─── Registration ─────────────────────────────────────────────────────────────

def register(
    he_path: Path,
    psr_path: Path,
    out_dir: Path,
    pixel_size: float,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    reg = WsiReg2D("HE_to_PSR", str(out_dir))

    # PSR first so it becomes the fixed/reference modality
    reg.add_modality(
        "PSR",
        str(psr_path),
        pixel_size=pixel_size,
        prepro_dict=PSR_PREPRO,
    )
    reg.add_modality(
        "HE",
        str(he_path),
        pixel_size=pixel_size,
        prepro_dict=HE_PREPRO,
    )

    # Moving = HE, Fixed = PSR.  Three-stage pipeline in one call.
    # wsireg chains the transforms automatically; each stage initialises
    # from the previous result, so the non-linear step only needs to
    # capture residual local deformation after the rigid+affine pass.
    reg.add_reg_path(
        "HE",
        "PSR",
        reg_params=["rigid", "affine", "nl"],
    )

    print("[1/2] Running registration  rigid → affine → non-linear …")
    reg.register_images()

    print("[2/2] Writing registered H&E as OME-TIFF pyramid …")
    reg.transform_images(file_writer="ome.tiff")

    print(f"\nDone. Output written to: {out_dir.resolve()}/")


# ─── QC overlay ───────────────────────────────────────────────────────────────

def _read_wsi_level(slide_path: Path, level: int) -> np.ndarray:
    """Read one pyramid level from a WSI (SVS/NDPI/TIFF) as an RGB ndarray."""
    slide = openslide.OpenSlide(str(slide_path))
    actual = min(level, slide.level_count - 1)
    w, h = slide.level_dimensions[actual]
    region = slide.read_region((0, 0), actual, (w, h)).convert("RGB")
    return np.array(region)


def _read_ometiff_level(tiff_path: Path, level: int) -> np.ndarray:
    """Read one level from an OME-TIFF pyramid produced by wsireg."""
    import tifffile
    with tifffile.TiffFile(str(tiff_path)) as tif:
        series = tif.series[0]
        actual = min(level, len(series.levels) - 1)
        arr = series.levels[actual].asarray()

    # Normalise axis order to HWC
    if arr.ndim == 2:
        arr = np.stack([arr, arr, arr], axis=-1)
    elif arr.ndim == 3 and arr.shape[0] in (1, 3, 4):
        arr = np.moveaxis(arr, 0, -1)

    return arr[..., :3].astype(np.uint8)


def _resize_to(arr: np.ndarray, h: int, w: int) -> np.ndarray:
    return np.array(Image.fromarray(arr).resize((w, h), Image.LANCZOS))


def _checkerboard(img_a: np.ndarray, img_b: np.ndarray, n: int = 8) -> np.ndarray:
    """Alternate n×n tiles between two same-shape RGB arrays."""
    h, w = img_a.shape[:2]
    th, tw = h // n, w // n
    board = img_a.copy()
    for r in range(n):
        for c in range(n):
            if (r + c) % 2 == 0:
                board[r * th:(r + 1) * th, c * tw:(c + 1) * tw] = \
                    img_b[r * th:(r + 1) * th, c * tw:(c + 1) * tw]
    return board


def qc_overlay(
    he_path: Path,
    psr_path: Path,
    reg_he_path: Path,
    out_dir: Path,
) -> None:
    """Save a four-panel QC PNG: PSR | H&E original | H&E registered | checkerboard."""
    print("Generating QC overlay …")

    psr_arr     = _read_wsi_level(psr_path, QC_PYRAMID_LEVEL)
    orig_he_arr = _read_wsi_level(he_path,  QC_PYRAMID_LEVEL)

    try:
        reg_he_arr = _read_ometiff_level(reg_he_path, QC_PYRAMID_LEVEL)
    except Exception as exc:
        print(f"Warning: could not read registered OME-TIFF ({exc}). "
              "Showing original H&E in that panel instead.")
        reg_he_arr = orig_he_arr

    # Normalise all panels to the PSR thumbnail dimensions
    h, w = psr_arr.shape[:2]
    orig_he_arr = _resize_to(orig_he_arr, h, w)
    reg_he_arr  = _resize_to(reg_he_arr,  h, w)

    board = _checkerboard(psr_arr, reg_he_arr)

    fig, axes = plt.subplots(1, 4, figsize=(24, 6))
    panels = [
        (psr_arr,     "PSR (fixed / reference)"),
        (orig_he_arr, "H&E (original, unregistered)"),
        (reg_he_arr,  "H&E (registered → PSR space)"),
        (board,       "Checkerboard  PSR | H&E-reg"),
    ]
    for ax, (img, title) in zip(axes, panels):
        ax.imshow(img)
        ax.set_title(title, fontsize=10)
        ax.axis("off")

    plt.tight_layout()
    out = out_dir / "qc_overlay.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    print(f"QC overlay saved: {out}")


# ─── CLI ──────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Register H&E → PSR mouse liver WSIs (rigid → affine → non-linear).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--he",          required=True,
                        help="Path to H&E WSI (.svs / .ndpi / .tiff / .scn)")
    parser.add_argument("--psr",         required=True,
                        help="Path to PSR WSI (.svs / .ndpi / .tiff / .scn)")
    parser.add_argument("--out",         default="./registration_output",
                        help="Output directory")
    parser.add_argument("--pixel-size",  type=float, default=DEFAULT_PIXEL_SIZE_UM,
                        help="Full-resolution pixel size in µm (check your scanner metadata)")
    parser.add_argument("--qc",          action="store_true",
                        help="Generate a QC checkerboard overlay PNG after registration")
    args = parser.parse_args()

    he_path  = Path(args.he)
    psr_path = Path(args.psr)
    out_dir  = Path(args.out)

    for p in (he_path, psr_path):
        if not p.exists():
            sys.exit(f"File not found: {p}")

    register(he_path, psr_path, out_dir, args.pixel_size)

    if args.qc:
        # wsireg names the output after the project + reg path
        candidates = sorted(out_dir.glob("*HE*registered*.ome.tiff"))
        if not candidates:
            candidates = sorted(out_dir.glob("*.ome.tiff"))
        if candidates:
            qc_overlay(he_path, psr_path, candidates[0], out_dir)
        else:
            print("Warning: could not find a registered OME-TIFF for QC — "
                  "check the output directory manually.")


if __name__ == "__main__":
    main()
