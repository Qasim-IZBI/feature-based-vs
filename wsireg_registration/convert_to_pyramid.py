#!/usr/bin/env python3
"""
Convert a flat .tif to a pyramidal OME-TIFF readable by wsireg.

Usage:
    python convert_to_pyramid.py input.tif
    python convert_to_pyramid.py input.tif --out output.ome.tiff
    python convert_to_pyramid.py input.tif --levels 5 --tile 256
"""

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import tifffile
except ImportError:
    sys.exit("tifffile not found. Install with: pip install tifffile")


def build_pyramid(img: np.ndarray, n_levels: int) -> list[np.ndarray]:
    """Downsample by 2× at each level using a simple average pooling."""
    levels = [img]
    for _ in range(n_levels - 1):
        prev = levels[-1]
        h, w = prev.shape[:2]
        # Trim to even dimensions before pooling
        prev = prev[: h - h % 2, : w - w % 2]
        if prev.ndim == 3:
            downsampled = (
                prev[0::2, 0::2].astype(np.uint32)
                + prev[1::2, 0::2]
                + prev[0::2, 1::2]
                + prev[1::2, 1::2]
            ) // 4
        else:
            downsampled = (
                prev[0::2, 0::2].astype(np.uint32)
                + prev[1::2, 0::2]
                + prev[0::2, 1::2]
                + prev[1::2, 1::2]
            ) // 4
        levels.append(downsampled.astype(img.dtype))
        if min(downsampled.shape[:2]) < 512:
            break
    return levels


def convert(src: Path, dst: Path, n_levels: int, tile_size: int) -> None:
    print(f"Reading  {src} …")
    img = tifffile.imread(str(src))
    print(f"Shape: {img.shape}  dtype: {img.dtype}")

    if img.ndim not in (2, 3):
        sys.exit(f"Unexpected image shape {img.shape}. Expected HW or HWC.")

    print(f"Building {n_levels}-level pyramid …")
    levels = build_pyramid(img, n_levels)
    for i, lvl in enumerate(levels):
        print(f"  level {i}: {lvl.shape}")

    print(f"Writing  {dst} …")
    with tifffile.TiffWriter(str(dst), bigtiff=True) as tw:
        options = dict(
            photometric="rgb" if img.ndim == 3 and img.shape[2] == 3 else "minisblack",
            tile=(tile_size, tile_size),
            compression="deflate",
            metadata=None,
        )
        # Full-resolution page with SubIFDs for pyramid levels below it
        tw.write(levels[0], subfiletype=0, **options,
                 subifds=len(levels) - 1)
        for lvl in levels[1:]:
            tw.write(lvl, subfiletype=1, **options)

    print(f"Done. Output: {dst}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert a flat .tif to a pyramidal OME-TIFF for wsireg.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("input",          help="Flat .tif input file")
    parser.add_argument("--out",          help="Output path (default: <input>.ome.tiff)")
    parser.add_argument("--levels", type=int, default=6,
                        help="Number of pyramid levels")
    parser.add_argument("--tile",   type=int, default=512,
                        help="Tile size in pixels (must be a multiple of 16)")
    args = parser.parse_args()

    src = Path(args.input)
    if not src.exists():
        sys.exit(f"File not found: {src}")

    dst = Path(args.out) if args.out else src.with_suffix(".ome.tiff")

    if args.tile % 16 != 0:
        sys.exit("--tile must be a multiple of 16 (TIFF spec requirement).")

    convert(src, dst, args.levels, args.tile)


if __name__ == "__main__":
    main()
