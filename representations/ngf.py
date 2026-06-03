"""
Normalized Gradient Field (NGF)
================================
Introduced by Haber & Modersitzki (2006) and adopted in HistokatFusion
(ANHIR winner) for cross-modal WSI registration.

Formula:
    NGF(I) = ∇I / (|∇I| + ε)

By dividing by the gradient magnitude, intensity information is discarded
and only gradient *orientation* is preserved. This makes the representation
robust to the large intensity differences between staining modalities
(H&E purple vs IHC brown vs PSR red), which is exactly the cross-modal
invariance property needed for conditioning virtual staining models.

Outputs:
    - NGF magnitude  : how strong the gradient is (edge strength map)
    - NGF orientation: gradient direction encoded as HSV colour wheel
    - NGF x / y      : raw normalised components (for use as channels)
"""

import numpy as np
from scipy.ndimage import sobel
from skimage.color import hsv2rgb


class NGF:
    """Compute the Normalized Gradient Field of a histology image."""

    def __init__(self, epsilon: float = 1e-4):
        """
        Parameters
        ----------
        epsilon : small constant to avoid division by zero (default 1e-4)
        """
        self.epsilon = epsilon

    def compute(self, rgb: np.ndarray) -> dict:
        """
        Parameters
        ----------
        rgb : np.ndarray  shape (H, W, 3)  dtype uint8

        Returns
        -------
        dict with keys:
            'gray'            – greyscale input  (H, W) float32 in [0,1]
            'gradient_mag'    – |∇I| before normalization, normalised to [0,1]
            'ngf_magnitude'   – |NGF| = |∇I|/(|∇I|+ε), values in [0,1)
            'ngf_x'           – x-component of NGF  in [-1, 1]
            'ngf_y'           – y-component of NGF  in [-1, 1]
            'ngf_orientation' – (H, W, 3) RGB, gradient direction as colour
            'ngf_magnitude_cm'– (H, W, 3) RGB, magnitude as colormap image
        """
        gray = self._to_gray(rgb)

        # Sobel gradients (scipy.ndimage.sobel uses a 3x3 kernel)
        gx = sobel(gray, axis=1).astype(np.float32)
        gy = sobel(gray, axis=0).astype(np.float32)

        mag = np.sqrt(gx ** 2 + gy ** 2)

        eps = self.epsilon
        ngf_x = gx / (mag + eps)
        ngf_y = gy / (mag + eps)
        ngf_mag = mag / (mag + eps)  # in [0, 1)

        # Orientation visualisation: angle → hue, magnitude → value
        angle = np.arctan2(ngf_y, ngf_x)           # [-π, π]
        hue = (angle + np.pi) / (2 * np.pi)        # [0, 1]
        sat = np.ones_like(hue)
        val = ngf_mag
        hsv = np.stack([hue, sat, val], axis=-1)
        orientation_rgb = hsv2rgb(hsv).astype(np.float32)

        # Magnitude colourmap (inferno)
        magnitude_cm = self._apply_colormap(ngf_mag, cmap="inferno")

        # Raw magnitude (for display without cmap)
        grad_mag_norm = mag / (mag.max() + 1e-8)

        return {
            "gray": gray,
            "gradient_mag": grad_mag_norm,
            "ngf_magnitude": ngf_mag,
            "ngf_x": ngf_x,
            "ngf_y": ngf_y,
            "ngf_orientation": orientation_rgb,
            "ngf_magnitude_cm": magnitude_cm,
        }

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _to_gray(rgb: np.ndarray) -> np.ndarray:
        """Convert uint8 RGB to float32 greyscale in [0, 1]."""
        r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
        gray = 0.2989 * r + 0.5870 * g + 0.1140 * b
        return (gray / 255.0).astype(np.float32)

    @staticmethod
    def _apply_colormap(x: np.ndarray, cmap: str = "inferno") -> np.ndarray:
        """Apply a matplotlib colormap and return an (H,W,3) float32 image."""
        import matplotlib.cm as mcm
        cm = mcm.get_cmap(cmap)
        return cm(x)[..., :3].astype(np.float32)
