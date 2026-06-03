"""
run_all.py — Cross-Modal Representation Inspector
===================================================
Generates all alternative image representations for one or more histology
images and saves comparison grids for visual inspection.

Usage
-----
Single image:
    python run_all.py --images path/to/HE.png --output ./results

Multiple images (cross-modal comparison):
    python run_all.py --images HE.png IHC.png PSR.png --output ./results

Force stain mode for color deconvolution:
    python run_all.py --images IHC.png --deconv hdab --output ./results

Output
------
For each image:
    results/<stem>_representations.png   – 2×4 grid of all representations

When multiple images are given:
    results/comparison_<rep>.png         – one image per representation,
                                           showing all input images side-by-side

Representations included
------------------------
1. Original (RGB)
2. Greyscale
3. NGF magnitude          (Normalized Gradient Field — cross-modal invariant)
4. NGF orientation        (colour-coded gradient direction)
5. DoG composite          (multi-scale Difference of Gaussians)
6. Phase Congruency       (intensity-invariant edge/feature map)
7. HOG fine (8×8)         (dense gradient orientation, fine scale)
8. HOG coarse (16×16)     (dense gradient orientation, coarse scale)
9. K-means K=3            (false-colour tissue compartment clustering)
10. K-means K=5           (finer compartment clustering)
11. Color Deconv ch0      (H&E: Haematoxylin | IHC: Haematoxylin)
12. Color Deconv ch1      (H&E: Eosin         | IHC: DAB marker)
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from PIL import Image

# Add parent directory to path so `representations` package is importable
# whether run from Feature_based/ or from representations/
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from representations.color_deconvolution import ColorDeconvolution
from representations.kmeans_color import KMeansColor
from representations.ngf import NGF
from representations.dog import DoG
from representations.phase_congruency import PhaseCongruency
from representations.hog_map import HOGMap


# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────

def load_image(path: str) -> np.ndarray:
    """Load an image and return as uint8 RGB numpy array."""
    img = Image.open(path).convert("RGB")
    return np.array(img, dtype=np.uint8)


def to_display(x) -> np.ndarray:
    """Convert any float/uint representation to uint8 RGB for display."""
    if x is None:
        return np.zeros((8, 8, 3), dtype=np.uint8)
    x = np.asarray(x, dtype=np.float32)
    if x.ndim == 2:                          # greyscale → RGB
        x = np.stack([x, x, x], axis=-1)
    x = np.clip(x, 0, 1)
    return (x * 255).astype(np.uint8)


def save_grid(panels: list, title: str, out_path: str, cols: int = 4):
    """
    Save a grid of labelled image panels.

    panels : list of (label, image_array_uint8_RGB)
    """
    rows = int(np.ceil(len(panels) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 4))
    axes = np.array(axes).flatten()

    for i, (ax, (label, img)) in enumerate(zip(axes, panels)):
        ax.imshow(img)
        ax.set_title(label, fontsize=9, pad=4)
        ax.axis("off")

    # Hide unused axes
    for ax in axes[len(panels):]:
        ax.axis("off")

    fig.suptitle(title, fontsize=13, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def save_comparison_grid(
    image_data: list,    # list of (stem, panels_dict)
    rep_key: str,
    rep_label: str,
    out_path: str,
):
    """One row per representation key, one column per image."""
    n = len(image_data)
    fig, axes = plt.subplots(1, n, figsize=(n * 4, 4))
    if n == 1:
        axes = [axes]

    for ax, (stem, panels_dict) in zip(axes, image_data):
        img = panels_dict.get(rep_key)
        if img is not None:
            ax.imshow(to_display(img))
        ax.set_title(stem, fontsize=9, pad=4)
        ax.axis("off")

    fig.suptitle(rep_label, fontsize=12, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


# ──────────────────────────────────────────────────────────────────────
# Per-image computation
# ──────────────────────────────────────────────────────────────────────

def compute_all_representations(rgb: np.ndarray, deconv_mode: str = "auto") -> dict:
    """Run all representations and return a flat dict of named arrays."""
    results = {}

    # 1. NGF
    ngf = NGF().compute(rgb)
    results["Original"]          = rgb.astype(np.float32) / 255.0
    results["Greyscale"]         = ngf["gray"]
    results["NGF Magnitude"]     = ngf["ngf_magnitude_cm"]
    results["NGF Orientation"]   = ngf["ngf_orientation"]

    # 2. DoG
    dog_result = DoG().compute(rgb)
    results["DoG Composite"]     = dog_result["composite_cm"]
    for label, d in zip(dog_result["scale_labels"], dog_result["dogs_norm"]):
        results[label]           = d

    # 3. Phase Congruency
    pc_result = PhaseCongruency().compute(rgb)
    results["Phase Congruency"]  = pc_result["pc_map_cm"]
    results["PC Binary Edges"]   = pc_result["pc_thresh"].astype(np.float32) / 255.0

    # 4. HOG
    hog_result = HOGMap().compute(rgb)
    results["HOG fine (8×8)"]   = hog_result["hog_fine_cm"]
    results["HOG coarse (16×16)"] = hog_result["hog_coarse_cm"]

    # 5. K-means
    km_result = KMeansColor(k_values=(3, 5)).compute(rgb)
    results["K-means K=3"]       = km_result[3]["false_colour"].astype(np.float32) / 255.0
    results["K-means K=5"]       = km_result[5]["false_colour"].astype(np.float32) / 255.0

    # 6. Color Deconvolution (H&E and/or H-DAB)
    cd = ColorDeconvolution(mode=deconv_mode)
    cd_result = cd.compute(rgb)
    mode_label = cd_result["mode"].upper()
    ch0_label  = cd_result.get("channel_0_label", "Ch0")
    ch1_label  = cd_result.get("channel_1_label", "Ch1")
    results[f"Deconv {ch0_label} ({mode_label})"] = cd_result["channel_0_rgb"]
    results[f"Deconv {ch1_label} ({mode_label})"] = cd_result["channel_1_rgb"]
    results[f"Deconv Overlay ({mode_label})"]     = cd_result["overlay"]

    # If auto mode detected H&E but user may also want H-DAB (for inspection),
    # also run H-DAB explicitly so both are always available
    if deconv_mode == "auto" and cd_result["mode"] == "he":
        cd_hdab = ColorDeconvolution(mode="hdab").compute(rgb)
        results["Deconv Haematoxylin (HDAB)"] = cd_hdab["channel_0_rgb"]
        results["Deconv DAB (HDAB)"]          = cd_hdab["channel_1_rgb"]
        results["Deconv Overlay (HDAB)"]      = cd_hdab["overlay"]

    return results


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Compute and visualise cross-modal image representations."
    )
    parser.add_argument(
        "--images", nargs="+", required=True,
        help="One or more image paths (H&E, IHC, PSR, etc.)"
    )
    parser.add_argument(
        "--output", default="./results",
        help="Output directory (created if it does not exist)"
    )
    parser.add_argument(
        "--deconv", default="auto",
        choices=["auto", "he", "hdab", "grayscale"],
        help="Color deconvolution mode (default: auto-detect)"
    )
    parser.add_argument(
        "--cols", type=int, default=4,
        help="Columns in per-image grid (default: 4)"
    )
    args = parser.parse_args()

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_image_data = []   # list of (stem, rep_dict) for cross-modal grid

    for img_path in args.images:
        stem = Path(img_path).stem
        print(f"\nProcessing: {img_path}")

        rgb = load_image(img_path)
        print(f"  Image size: {rgb.shape[1]}×{rgb.shape[0]} px")

        reps = compute_all_representations(rgb, deconv_mode=args.deconv)
        all_image_data.append((stem, reps))

        # Build panel list for per-image grid
        panels = [(k, to_display(v)) for k, v in reps.items()]

        out_file = str(out_dir / f"{stem}_representations.png")
        save_grid(panels, title=f"Representations — {stem}",
                  out_path=out_file, cols=args.cols)

    # ── Cross-image comparison grids (one file per representation) ──
    if len(all_image_data) > 1:
        print("\nGenerating cross-image comparison grids …")
        # Collect all representation keys that appear in every image
        all_keys = list(all_image_data[0][1].keys())
        comp_dir = out_dir / "comparisons"
        comp_dir.mkdir(exist_ok=True)

        for key in all_keys:
            safe_key = key.replace(" ", "_").replace("/", "-").replace("×", "x")
            out_file = str(comp_dir / f"compare_{safe_key}.png")
            save_comparison_grid(
                image_data=all_image_data,
                rep_key=key,
                rep_label=key,
                out_path=out_file,
            )

    print(f"\nDone. Results saved to: {out_dir}")


if __name__ == "__main__":
    main()
