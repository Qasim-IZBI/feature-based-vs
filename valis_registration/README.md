# H&E ↔ PSR WSI Registration — VALIS

Registers H&E and PSR mouse liver WSIs using VALIS (feature-based, SIFT).

## Why VALIS over wsireg/Elastix

| | wsireg (Elastix) | VALIS |
|---|---|---|
| Method | Intensity-based (mutual info) | Feature-based (SIFT) |
| Image size handling | Registers at full/high res → slow | Always uses thumbnail (`MAX_PROC_SIZE`) → fast |
| Multi-stain robustness | Needs careful channel tuning | SIFT works across stain types natively |
| Flat .tif support | Needs pyramid conversion first | Reads flat .tif directly |
| Speed (20k×20k pair) | 100+ min (hangs on output) | 3–10 min |

## Setup

### 1. Install system libraries

```bash
# macOS
brew install vips openslide

# Ubuntu / Debian
sudo apt install libvips-dev libopenslide-dev
```

### 2. Install Python dependencies

```bash
pip install -r requirements.txt
```

## Usage

Open `registration_pipeline.ipynb`, set `HE_DIR`, `PSR_DIR`, and `OUTPUT_DIR`
in the config cell, then run all cells.

Key settings:

| Parameter | Default | Effect |
|---|---|---|
| `MAX_PROC_SIZE` | `850` | Thumbnail size for registration — increase for more accuracy, decrease for speed |
| `NON_RIGID` | `False` | Add optical-flow non-rigid step after rigid+affine — try only if QC overlay shows local warping |

## Output layout

```
valis_output/
├── mouse01/
│   ├── registered/
│   │   ├── mouse01_HE.ome.tiff     ← H&E (reference, unchanged)
│   │   └── mouse01_PSR.ome.tiff    ← PSR registered into H&E space
│   └── ...                         ← VALIS intermediate files
└── qc/
    └── mouse01_qc.png              ← four-panel QC overlay
```
