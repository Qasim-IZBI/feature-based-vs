"""
Phase Congruency
================
Kovesi (1999) implementation using log-Gabor filters in the frequency domain.

Phase congruency detects features at locations where the Fourier components
of the image are maximally in phase. Unlike gradient-based methods, it is
entirely *intensity-invariant*: it does not depend on the absolute brightness
or contrast of the image — only on the local phase structure.

This makes it ideal for cross-modal histology comparison:
  - A gland boundary in H&E and the same boundary in IHC/PSR will both
    produce a strong phase congruency response, even though the intensities
    are completely different.

Reference:
  Kovesi, P. (1999). Image features from phase congruency.
  Videre: Journal of Computer Vision Research, 1(3).

Parameters (defaults follow Kovesi's recommendations):
  nscale       : number of filter scales              (default 4)
  norient      : number of filter orientations         (default 6)
  minWaveLength: wavelength of smallest filter (px)   (default 3)
  mult         : scaling factor between successive filters (default 2.1)
  sigmaOnf     : bandwidth of each filter (log-Gabor) (default 0.55)
  k            : noise compensation factor            (default 2.0)
  cutOff       : low-pass weighting cut-off           (default 0.5)
  g            : sharpness of cut-off                 (default 10)
"""

import numpy as np


class PhaseCongruency:
    """Compute phase congruency of a greyscale or RGB histology image."""

    def __init__(
        self,
        nscale: int = 4,
        norient: int = 6,
        min_wave_length: int = 3,
        mult: float = 2.1,
        sigma_onf: float = 0.55,
        k: float = 2.0,
        cut_off: float = 0.5,
        g: int = 10,
    ):
        self.nscale = nscale
        self.norient = norient
        self.min_wave_length = min_wave_length
        self.mult = mult
        self.sigma_onf = sigma_onf
        self.k = k
        self.cut_off = cut_off
        self.g = g

    def compute(self, rgb: np.ndarray) -> dict:
        """
        Parameters
        ----------
        rgb : np.ndarray  shape (H, W, 3)  dtype uint8

        Returns
        -------
        dict with keys:
            'pc_map'      – (H, W) float32, phase congruency values in [0,1]
            'pc_map_cm'   – (H, W, 3) RGB with inferno colormap
            'orientation' – (H, W) float32, dominant orientation in [0, π]
            'pc_thresh'   – (H, W) uint8 binary edge map (Otsu threshold)
        """
        gray = self._to_gray(rgb)
        pc, orientation = self._phase_congruency(gray)
        pc_cm = self._apply_colormap(pc, "inferno")

        # Binary edge map via Otsu threshold
        from skimage.filters import threshold_otsu
        thresh = threshold_otsu(pc)
        pc_thresh = (pc > thresh).astype(np.uint8) * 255

        return {
            "pc_map": pc,
            "pc_map_cm": pc_cm,
            "orientation": orientation,
            "pc_thresh": pc_thresh,
        }

    # ------------------------------------------------------------------
    # Core algorithm (Kovesi 1999)
    # ------------------------------------------------------------------

    def _phase_congruency(self, gray: np.ndarray):
        """
        Compute phase congruency from a greyscale float32 image.
        Returns (pc_map, orientation_map).
        """
        H, W = gray.shape
        eps = 1e-10

        # ----- build frequency-domain grids -----
        # u, v in [-0.5, 0.5]
        u = np.fft.fftfreq(W).astype(np.float32)
        v = np.fft.fftfreq(H).astype(np.float32)
        U, V = np.meshgrid(u, v)
        radius = np.sqrt(U ** 2 + V ** 2)
        radius[0, 0] = 1.0   # avoid log(0)
        theta_grid = np.arctan2(V, U)  # orientation of each freq component

        # Pre-compute log-Gabor radial component at each scale
        log_gabor_scales = []
        wave_length = self.min_wave_length
        for _ in range(self.nscale):
            fo = 1.0 / wave_length
            log_gabor = np.exp(
                -(np.log(radius / fo)) ** 2 / (2 * np.log(self.sigma_onf) ** 2)
            )
            log_gabor[0, 0] = 0.0   # remove DC
            log_gabor_scales.append(log_gabor)
            wave_length *= self.mult

        # FFT of input
        IM = np.fft.fft2(gray.astype(np.float64))

        # Accumulators across orientations
        total_sum_an = np.zeros((H, W), dtype=np.float64)
        total_energy = np.zeros((H, W), dtype=np.float64)
        total_mean_psi = np.zeros((H, W), dtype=np.complex128)

        for o in range(self.norient):
            angle_o = o * np.pi / self.norient   # orientation of this filter
            # Angular spread (half-width)
            theta_sigma = np.pi / (self.norient * 2)

            # Compute angular difference (wrap to [-π/2, π/2])
            dtheta = theta_grid - angle_o
            dtheta = np.abs(np.arctan2(np.sin(dtheta), np.cos(dtheta)))

            # Angular Gaussian component
            spread = np.exp(-dtheta ** 2 / (2 * theta_sigma ** 2))

            # Accumulators for this orientation
            sum_an = np.zeros((H, W), dtype=np.float64)
            sum_an_cos = np.zeros((H, W), dtype=np.float64)
            sum_an_sin = np.zeros((H, W), dtype=np.float64)
            sum_e = np.zeros((H, W), dtype=np.float64)
            sum_o = np.zeros((H, W), dtype=np.float64)

            for s in range(self.nscale):
                filt = log_gabor_scales[s] * spread
                # Apply filter
                EO = np.fft.ifft2(IM * filt)
                An = np.abs(EO)
                phase = np.angle(EO)

                sum_an += An
                sum_an_cos += An * np.cos(phase)
                sum_an_sin += An * np.sin(phase)
                sum_e += EO.real
                sum_o += EO.imag

            # Mean phase angle
            mean_cos = sum_an_cos / (sum_an + eps)
            mean_sin = sum_an_sin / (sum_an + eps)
            mean_psi = np.arctan2(mean_sin, mean_cos)

            # Phase deviation from mean
            for s in range(self.nscale):
                filt = log_gabor_scales[s] * spread
                EO = np.fft.ifft2(IM * filt)
                An = np.abs(EO)
                phase = np.angle(EO)
                delta_phi = phase - mean_psi
                # Weighted phase: cos(delta_phi) - |sin(delta_phi)|
                w = An
                total_energy += w * (np.cos(delta_phi) - np.abs(np.sin(delta_phi)))
                total_sum_an += An

            total_mean_psi += np.exp(1j * mean_psi)

        # ----- noise compensation -----
        # Median estimate of noise energy (from the smallest scale)
        noise_power = np.median(np.abs(
            np.fft.ifft2(IM * log_gabor_scales[0])
        ))
        tau = noise_power * np.sqrt(np.pi / 2)
        noise_threshold = (
            tau * self.k * np.sqrt(2 * np.log(1.0 / (1 - self.cut_off ** self.g)))
            if self.cut_off > 0 else 0.0
        )
        # Soft threshold on energy
        energy_thresh = np.maximum(total_energy - noise_threshold, 0.0)

        # Low-pass weighting (suppress very low-energy regions)
        weight = 1.0 / (1.0 + np.exp(
            self.g * (self.cut_off - total_sum_an / (self.nscale * self.norient + eps))
        ))

        pc = weight * energy_thresh / (total_sum_an + eps)
        pc = np.clip(pc, 0, 1).astype(np.float32)

        # Dominant orientation
        orientation = np.angle(total_mean_psi).astype(np.float32)  # [-π, π]
        orientation = (orientation % np.pi).astype(np.float32)      # [0, π]

        return pc, orientation

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _to_gray(rgb: np.ndarray) -> np.ndarray:
        r, g, b = rgb[..., 0], rgb[..., 1], rgb[..., 2]
        gray = 0.2989 * r + 0.5870 * g + 0.1140 * b
        return (gray / 255.0).astype(np.float32)

    @staticmethod
    def _apply_colormap(x: np.ndarray, cmap: str = "inferno") -> np.ndarray:
        import matplotlib.cm as mcm
        cm = mcm.get_cmap(cmap)
        return cm(x)[..., :3].astype(np.float32)
