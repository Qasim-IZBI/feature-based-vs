"""
Color Deconvolution
===================
Separates histological images into individual stain channels using the
Ruifrok & Johnston (2001) optical density matrix method.

Supported stain systems:
  - H&E  : Haematoxylin (nuclei, blue-purple) + Eosin (cytoplasm/ECM, pink)
  - H-DAB: Haematoxylin (nuclei, blue) + DAB (protein marker, brown)
            This is the standard IHC stain system.

For stains that do not fit either system (e.g. PSR, PAS), a fallback
grayscale + inverted-grayscale pair is returned (STAR-style normalization).
"""

import numpy as np
from skimage.color import separate_stains, hed_from_rgb


# Stain vectors from Ruifrok & Johnston (2001) and scikit-image defaults.
# Each row is [R, R, B] optical density contribution of one stain.
#
# hed_from_rgb covers all three:  H=row0, E=row1, DAB=row2
#
# For H-DAB we use H (row0) and DAB (row2) and discard E (row1).

# Convenience alias — scikit-image ships this matrix
HED_MATRIX = hed_from_rgb  # shape (3,3)


class ColorDeconvolution:
    """Separate a stained image into its constituent stain channels."""

    SUPPORTED_MODES = ("he", "hdab", "auto", "grayscale")

    def __init__(self, mode: str = "auto"):
        """
        Parameters
        ----------
        mode : str
            'he'        – force H&E deconvolution
            'hdab'      – force H-DAB deconvolution (IHC)
            'auto'      – guess from dominant stain colour (default)
            'grayscale' – skip deconvolution, return grayscale + inverted
        """
        assert mode in self.SUPPORTED_MODES, (
            f"mode must be one of {self.SUPPORTED_MODES}, got '{mode}'"
        )
        self.mode = mode

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def compute(self, rgb: np.ndarray) -> dict:
        """
        Parameters
        ----------
        rgb : np.ndarray  shape (H, W, 3)  dtype uint8

        Returns
        -------
        dict with keys:
            'mode'          – which deconvolution was applied
            'channel_0'     – first stain channel  (H&E: H | H-DAB: H | gray: gray)
            'channel_1'     – second stain channel (H&E: E | H-DAB: DAB | gray: inv)
            'channel_2'     – third channel if available (H&E only: residual)
            'channel_0_rgb' – channel_0 reconstructed in pseudo-colour
            'channel_1_rgb' – channel_1 reconstructed in pseudo-colour
            'overlay'       – false-colour overlay of ch0 (red) + ch1 (cyan)
        """
        mode = self._resolve_mode(rgb) if self.mode == "auto" else self.mode

        if mode in ("he", "hdab"):
            return self._deconvolve(rgb, mode)
        else:
            return self._grayscale_fallback(rgb)

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _resolve_mode(rgb: np.ndarray) -> str:
        """Guess stain mode from the dominant hue of the image."""
        # Convert to float and look at mean hue in HSV space
        from skimage.color import rgb2hsv
        hsv = rgb2hsv(rgb.astype(np.float32) / 255.0)
        hue = hsv[..., 0]
        sat = hsv[..., 1]
        # Only consider saturated pixels (tissue, not white background)
        tissue_mask = sat > 0.1
        if tissue_mask.sum() == 0:
            return "grayscale"
        mean_hue = hue[tissue_mask].mean()  # in [0, 1]
        mean_hue_deg = mean_hue * 360.0

        # H&E: mix of purple (270°) and pink (330–350°) → mean ~300°
        # H-DAB (IHC): mix of blue (240°) and brown (30°) → mean varies
        # PSR / unknown: anything else
        if 250 < mean_hue_deg < 360 or mean_hue_deg < 30:
            # Purple/pink/red dominant → likely H&E or IHC
            # Distinguish by presence of brown (DAB) vs pink (eosin)
            # DAB is in hue range 20–50°
            brown_frac = ((hue[tissue_mask] * 360 > 15) &
                          (hue[tissue_mask] * 360 < 55)).mean()
            if brown_frac > 0.05:
                return "hdab"
            return "he"
        return "grayscale"

    @staticmethod
    def _deconvolve(rgb: np.ndarray, mode: str) -> dict:
        img = rgb.astype(np.float32) / 255.0
        # separate_stains returns optical density concentrations [H, W, 3]
        stains = separate_stains(img, HED_MATRIX)  # channels: H, E, DAB

        if mode == "he":
            ch0 = stains[..., 0]  # Haematoxylin
            ch1 = stains[..., 1]  # Eosin
            ch2 = stains[..., 2]  # Residual / DAB
            ch0_label, ch1_label = "Haematoxylin", "Eosin"
        else:  # hdab
            ch0 = stains[..., 0]  # Haematoxylin
            ch1 = stains[..., 2]  # DAB
            ch2 = stains[..., 1]  # Eosin (residual for IHC)
            ch0_label, ch1_label = "Haematoxylin", "DAB"

        def _norm(x):
            """Normalize to [0, 1] for display."""
            x = x - x.min()
            r = x.max()
            return x / r if r > 1e-8 else x

        ch0_n = _norm(ch0)
        ch1_n = _norm(ch1)

        # Pseudo-colour: stain 0 → purple tint, stain 1 → pink/brown tint
        ch0_rgb = ColorDeconvolution._channel_to_rgb(ch0_n, colour=(0.55, 0.20, 0.80))
        ch1_rgb = ColorDeconvolution._channel_to_rgb(ch1_n, colour=(0.90, 0.20, 0.20))

        # Simple overlay: ch0 in red channel, ch1 in green channel
        overlay = np.stack([ch0_n, ch1_n, np.zeros_like(ch0_n)], axis=-1)
        overlay = np.clip(overlay, 0, 1)

        return {
            "mode": mode,
            "channel_0": ch0_n,
            "channel_1": ch1_n,
            "channel_2": _norm(ch2),
            "channel_0_rgb": ch0_rgb,
            "channel_1_rgb": ch1_rgb,
            "overlay": overlay,
            "channel_0_label": ch0_label,
            "channel_1_label": ch1_label,
        }

    @staticmethod
    def _grayscale_fallback(rgb: np.ndarray) -> dict:
        from skimage.color import rgb2gray
        gray = rgb2gray(rgb.astype(np.float32) / 255.0)
        inv = 1.0 - gray
        return {
            "mode": "grayscale",
            "channel_0": gray,
            "channel_1": inv,
            "channel_2": None,
            "channel_0_rgb": np.stack([gray, gray, gray], axis=-1),
            "channel_1_rgb": np.stack([inv, inv, inv], axis=-1),
            "overlay": np.stack([inv, gray, np.zeros_like(gray)], axis=-1),
            "channel_0_label": "Grayscale",
            "channel_1_label": "Inverted",
        }

    @staticmethod
    def _channel_to_rgb(
        channel: np.ndarray, colour: tuple = (1.0, 0.0, 0.0)
    ) -> np.ndarray:
        """Map a single-channel [0,1] map to a tinted RGB image on white background."""
        r, g, b = colour
        rgb = np.ones((*channel.shape, 3), dtype=np.float32)
        rgb[..., 0] = 1.0 - channel * (1.0 - r)
        rgb[..., 1] = 1.0 - channel * (1.0 - g)
        rgb[..., 2] = 1.0 - channel * (1.0 - b)
        return np.clip(rgb, 0, 1)
