# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Research Context

This repository is an active research project exploring **registration-informed conditioning for virtual histology staining diffusion models**. The central question is whether classical, cross-modal-invariant image representations (derived from the WSI registration literature) can replace or improve the naive gradient-map conditioning used in MIU-Diff.

See `research_discussion.md` for the full discussion, motivation, and open questions. Read it before making any architectural decisions.

Key papers:
- `MIU-Diff/2506.23184v1.pdf` — the target diffusion model being improved
- `registration_papers/` — six WSI registration methods whose intermediate representations are being repurposed as conditioning signals

## Running the Representation Inspector

```bash
# Single image — generates a grid of all representations
python representations/run_all.py --images path/to/HE.png --output ./results

# Multi-image cross-modal comparison (most useful for research)
python representations/run_all.py --images HE.png IHC.png PSR.png --output ./results

# Force a specific color deconvolution mode
python representations/run_all.py --images IHC.png --deconv hdab --output ./results
# --deconv options: auto (default) | he | hdab | grayscale

# Control grid layout
python representations/run_all.py --images HE.png --cols 6 --output ./results
```

The script can be invoked from any directory; it auto-configures `sys.path`.

Dependencies: `numpy scipy scikit-image scikit-learn matplotlib Pillow`

## Architecture of `representations/`

Each module is a self-contained class with a single `.compute(rgb: np.ndarray) -> dict` method. Input is always `uint8 RGB (H, W, 3)`. All outputs are normalised to `[0, 1] float32` (scalar maps) or `uint8` (false-colour maps).

| Module | Class | What it computes | Cross-modal invariant? |
|---|---|---|---|
| `color_deconvolution.py` | `ColorDeconvolution` | H&E → {Haematoxylin, Eosin}; IHC → {Haematoxylin, DAB} | Partial (stain-specific) |
| `kmeans_color.py` | `KMeansColor` | Pixel clustering in LAB space at K=3 and K=5 | Partial |
| `ngf.py` | `NGF` | Normalized Gradient Field: `∇I / (|∇I| + ε)` | **Yes** |
| `dog.py` | `DoG` | Multi-scale Difference of Gaussians at σ pairs (1,2), (2,4), (4,8) | **Yes** |
| `phase_congruency.py` | `PhaseCongruency` | Kovesi (1999) log-Gabor phase congruency | **Yes** (intensity-invariant) |
| `hog_map.py` | `HOGMap` | Dense HOG at 8×8 (fine) and 16×16 (coarse) cell sizes | Partial |

`run_all.py` orchestrates all modules and produces:
- Per-image grids: `results/<stem>_representations.png`
- Cross-image comparison grids: `results/comparisons/compare_<rep>.png`

## Key Design Decisions

**Color deconvolution mode detection** (`color_deconvolution.py:_resolve_mode`): Auto-detection uses mean hue of tissue pixels (saturated pixels only, sat > 0.1). Hue 250–360°/0–30° → H&E or IHC; brown fraction (hue 15–55°) > 5% → H-DAB; anything else → grayscale fallback. **PSR, PAS, and unknown stains always fall back to grayscale+inverted** — there is no PSR-specific deconvolution. When auto-detects H&E, it *also* runs H-DAB automatically, because both are always useful for inspection. Override with `--deconv` if misdetected.

**NGF is the primary cross-modal signal of interest** (`ngf.py`). It is the direct replacement candidate for the Sobel gradient map in MIU-Diff's MI estimator. The key property: only orientation is preserved, intensity is removed. For MIU-Diff integration, use the raw `ngf_x` and `ngf_y` output keys (normalized gradient components, range [-1, 1]) — not the display variants (`ngf_magnitude_cm`, `ngf_orientation`) which are colormap images for visual inspection only.

**Phase congruency is the most computationally expensive** (`phase_congruency.py`). It runs FFT-based log-Gabor filtering at `nscale × norient` = 24 filter banks by default. For large WSI patches, reduce `nscale` or `norient` first.

**Stain-specific conditioning strategy** (from `research_discussion.md`):
- H&E → IHC: nuclei are shared → Haematoxylin channel + NGF
- H&E → PSR: collagen architecture is shared, nuclei are NOT visible in PSR → tissue compartment maps + NGF
- For any unknown target stain: NGF and phase congruency are always safe choices

## Adding a New Representation

1. Create `representations/<name>.py` with a class exposing `.compute(rgb) -> dict`
2. Add import to `representations/__init__.py`
3. Call it inside `compute_all_representations()` in `run_all.py` and add results to the `results` dict
