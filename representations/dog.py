"""
Difference of Gaussians (DoG)
==============================
Computes the DoG response at multiple scale pairs. DoG is the underlying
detection mechanism of SIFT and approximates the Laplacian of Gaussian (LoG).

It highlights structurally stable locations across scale:
  - Gland boundaries
  - Vessel walls
  - Nuclear cluster borders
  - Tissue interface lines

Crucially, DoG responds to *structural transitions* independently of
absolute staining intensity, making it relatively modality-invariant
compared to raw gradient maps.

Scale pairs used by default: (1,2), (2,4), (4,8) pixels sigma.
The composite is the per-pixel maximum absolute response across all scales.
"""

import numpy as np
from scipy.ndimage import gaussian_filter


class DoG:
    """Multi-scale Difference of Gaussians on greyscale histology images."""

    # Default sigma pairs (fine → coarse)
    DEFAULT_SIGMA_PAIRS = [(1, 2), (2, 4), (4, 8)]

    def __init__(self, sigma_pairs: list = None):
        """
        Parameters
        ----------
        sigma_pairs : list of (sigma_fine, sigma_coarse) tuples
                      If None, uses DEFAULT_SIGMA_PAIRS.
        """
        self.sigma_pairs = sigma_pairs or self.DEFAULT_SIGMA_PAIRS

    def compute(self, rgb: np.ndarray) -> dict:
        """
        Parameters
        ----------
        rgb : np.ndarray  shape (H, W, 3)  dtype uint8

        Returns
        -------
        dict with keys:
            'gray'           – greyscale input (H, W) float32 in [0,1]
            'dogs'           – list of (H,W) float32 DoG maps, one per scale pair
            'dogs_norm'      – same, normalised to [0,1] for display
            'composite'      – (H,W) float32, max |DoG| across scales, norm [0,1]
            'composite_cm'   – (H,W,3) RGB composite with inferno colormap
            'scale_labels'   – list of strings, e.g. "DoG σ(1→2)"
        """
        gray = self._to_gray(rgb)

        dogs = []
        dogs_norm = []
        labels = []

        for s1, s2 in self.sigma_pairs:
            d = gaussian_filter(gray, s1) - gaussian_filter(gray, s2)
            dogs.append(d)
            dogs_norm.append(self._norm(np.abs(d)))
            labels.append(f"DoG σ({s1}→{s2})")

        # Composite: maximum absolute response across scales
        stack = np.stack([np.abs(d) for d in dogs], axis=0)
        composite = stack.max(axis=0)
        composite_norm = self._norm(composite)

        composite_cm = self._apply_colormap(composite_norm, cmap="inferno")

        return {
            "gray": gray,
            "dogs": dogs,
            "dogs_norm": dogs_norm,
            "composite": composite_norm,
            "composite_cm": composite_cm,
            "scale_labels": labels,
        }

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _to_gray(rgb: np.ndarray) -> np.ndarray:
        r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
        gray = 0.2989 * r + 0.5870 * g + 0.1140 * b
        return (gray / 255.0).astype(np.float32)

    @staticmethod
    def _norm(x: np.ndarray) -> np.ndarray:
        mn, mx = x.min(), x.max()
        if mx - mn < 1e-8:
            return np.zeros_like(x)
        return ((x - mn) / (mx - mn)).astype(np.float32)

    @staticmethod
    def _apply_colormap(x: np.ndarray, cmap: str = "inferno") -> np.ndarray:
        import matplotlib.cm as mcm
        cm = mcm.get_cmap(cmap)
        return cm(x)[..., :3].astype(np.float32)
