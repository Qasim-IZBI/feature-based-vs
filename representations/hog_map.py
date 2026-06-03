"""
Dense HOG Map
=============
Computes a dense Histogram of Oriented Gradients (HOG) feature map and
returns it as a visualizable image.

HOG describes local gradient orientation distributions rather than absolute
intensities, making it relatively insensitive to staining intensity differences
across modalities. It captures tissue texture and boundary structure in a
representation that can be compared across H&E, IHC, and PSR images.

Two cell sizes are computed:
  - Fine  (8×8  px): captures individual nucleus-level texture
  - Coarse (16×16 px): captures gland/tissue-region level structure
"""

import numpy as np
from skimage.feature import hog
from skimage.color import rgb2gray


class HOGMap:
    """Compute and visualise dense HOG for cross-modal comparison."""

    def __init__(
        self,
        pixels_per_cell_fine: tuple = (8, 8),
        pixels_per_cell_coarse: tuple = (16, 16),
        cells_per_block: tuple = (2, 2),
        orientations: int = 9,
    ):
        self.pixels_per_cell_fine = pixels_per_cell_fine
        self.pixels_per_cell_coarse = pixels_per_cell_coarse
        self.cells_per_block = cells_per_block
        self.orientations = orientations

    def compute(self, rgb: np.ndarray) -> dict:
        """
        Parameters
        ----------
        rgb : np.ndarray  shape (H, W, 3)  dtype uint8

        Returns
        -------
        dict with keys:
            'gray'          – greyscale input (H, W) float32 in [0,1]
            'hog_fine'      – (H, W) float32, HOG visualisation (8×8 cells)
            'hog_coarse'    – (H, W) float32, HOG visualisation (16×16 cells)
            'hog_fine_cm'   – (H, W, 3) RGB with viridis colormap
            'hog_coarse_cm' – (H, W, 3) RGB with viridis colormap
        """
        gray = rgb2gray(rgb.astype(np.float32) / 255.0)

        hog_fine = self._compute_hog(gray, self.pixels_per_cell_fine)
        hog_coarse = self._compute_hog(gray, self.pixels_per_cell_coarse)

        return {
            "gray": gray,
            "hog_fine": hog_fine,
            "hog_coarse": hog_coarse,
            "hog_fine_cm": self._apply_colormap(hog_fine, "viridis"),
            "hog_coarse_cm": self._apply_colormap(hog_coarse, "viridis"),
        }

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _compute_hog(self, gray: np.ndarray, pixels_per_cell: tuple) -> np.ndarray:
        _, hog_image = hog(
            gray,
            orientations=self.orientations,
            pixels_per_cell=pixels_per_cell,
            cells_per_block=self.cells_per_block,
            visualize=True,
            channel_axis=None,
        )
        # Normalise to [0, 1] for display
        mn, mx = hog_image.min(), hog_image.max()
        if mx - mn > 1e-8:
            hog_image = (hog_image - mn) / (mx - mn)
        return hog_image.astype(np.float32)

    @staticmethod
    def _apply_colormap(x: np.ndarray, cmap: str = "viridis") -> np.ndarray:
        import matplotlib.cm as mcm
        cm = mcm.get_cmap(cmap)
        return cm(x)[..., :3].astype(np.float32)
